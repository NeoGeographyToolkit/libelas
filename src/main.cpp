// Modifications by Oleg Alexandrov @ NASA Ames.
// See README.TXT for what changed. Released under the same license.

/*
Copyright 2011. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

This file is part of libelas.
Authors: Andreas Geiger

libelas is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or any later version.

libelas is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
libelas; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301, USA 
*/

#include <iostream>
#include <limits>
#include <algorithm>
#include <tiffio.h>
#include <cmath>
#include "elas.h"
#include "image.h"

extern "C" {
#include "iio.h"
}

using namespace std;

// compute disparities of pgm image input pair file_1, file_2
void process (const char* file_1, const char* file_2) {

  cout << "Processing: " << file_1 << ", " << file_2 << endl;

  int width, height; 
  float * limg = iio_read_image_float(file_1, &width, &height);

  int rwidth, rheight; 
  float * rimg = iio_read_image_float(file_2, &rwidth, &rheight);

  std::cout << "--deal with padding!" << std::endl;
  int pad = 45;

  // check for correct size
  if (width <= 0 || height <= 0 || rwidth <= 0 || rheight <= 0 ||
      width != rwidth || height != rheight) {
    std::cout << "ERROR: Images must be of same size, but got: " << std::endl;
    std::cout << "  left_image:  " << width <<  " x " << height  << std::endl;
    std::cout << "  right_image: " << rwidth <<  " x " << rheight << std::endl;
    return; 
  }

  int pad_width = width + pad;
  
  // Convert the image pixels from floats in [0, 1] to uint8_t in [0, 255].
  // To ensure a positive disparity, pad the left image on the left with enough zeros,
  // and pad the right image on the right with the same amount.
  uint8_t* l_img_pad = (uint8_t*)malloc(pad_width * height * sizeof(uint8_t));
  uint8_t* r_img_pad = (uint8_t*)malloc(pad_width * height * sizeof(uint8_t));

  // Process the left image
  for (int i = 0; i < pad_width * height; i++) 
    l_img_pad[i] = 0;

  int count = 0;
  for (int ih = 0; ih < height; ih++) {
    for (int iw = 0; iw < width; iw++) {

      float val = limg[count];
      if (std::isnan(val)) 
        val = 0.0;

      val *= 255.0;
      val = round(val);
      if (val < 0) 
        val = 0;
      if (val > 255.0)
        val = 255.0;

      // The zero padding is on the left, so add the 'pad' value
      l_img_pad[ih * pad_width + iw + pad] = (uint8_t)val;
      
      count++;
    }
  }

#if 0
  // For debugging
   char l_img_pad_file[] = "l_img_pad.tif";
   int ch = 1;
   iio_save_image_uint8_vec((char*)l_img_pad_file, l_img_pad, pad_width, height, ch);
#endif
   
  // Process the right image
  for (int i = 0; i < pad_width * height; i++) 
    r_img_pad[i] = 0;

  count = 0;
  for (int ih = 0; ih < height; ih++) {
    for (int iw = 0; iw < width; iw++) {

      float val = rimg[count];
      if (std::isnan(val)) 
        val = 0.0;

      val *= 255.0;
      val = round(val);
      if (val < 0) 
        val = 0;
      if (val > 255.0)
        val = 255.0;

      // The zero padding is on the right, so no need to add anything
      r_img_pad[ih * pad_width + iw] = (uint8_t)val;
      
      count++;
    }
  }
  
  // allocate memory for disparity images
  const int32_t dims[3] = {pad_width, height, pad_width}; // bytes per line = width
  float* lr_disp_pad = (float*)malloc(pad_width * height * sizeof(float));
  float* rl_disp_pad = (float*)malloc(pad_width * height * sizeof(float));
  
  // process
  Elas::parameters param;
  param.postprocess_only_left = false;
  Elas elas(param);
  elas.process(l_img_pad, r_img_pad, lr_disp_pad, rl_disp_pad, dims);

  // When the disparities have a big jump, the ones after the jump are outliers.
  // TODO(oalexan1): This is fragile.
  float* sorted_disp = (float*)malloc(pad_width * height * sizeof(float));
  for (int i = 0; i < pad_width * height; i++)
    sorted_disp[i] = lr_disp_pad[i];
  std::sort(sorted_disp, sorted_disp + pad_width * height);

  float max_jump = -1.0, max_valid = -1.0;
  for (int i = 0; i < pad_width * height - 1; i++) {
    float a = sorted_disp[i];
    float b = sorted_disp[i + 1];

    if (a < 0) 
      continue; // negative disparities are outliers

    if (b - a > max_jump) {
      max_jump = b - a;
      max_valid = a;
    }
  }

  free(sorted_disp);

  for (int i = 0; i < pad_width * height; i++)
    if (lr_disp_pad[i] > max_valid)
      lr_disp_pad[i] = -10.0; // invalidate the outlier disparity
  
  // Remove the padding, and subtract the padding value from the disparity
  // To undo the previous operations

  count = 0;
  float* lr_disp = (float*)malloc(width * height * sizeof(float));
  for (int ih = 0; ih < height; ih++) {
    for (int iw = 0; iw < width; iw++) {
      
      lr_disp[count] = lr_disp_pad[ih * pad_width + iw + pad];
      count++;
      
    }
  }
  
  // Negative disparities are set to NaN.
  // Flip the sign of the lr disparity as the libelas convention is
  // the opposite we expect.  Ignore the rl disparity. 
  float nan = std::numeric_limits<float>::quiet_NaN();
  for (int32_t i = 0; i < width * height; i++) {
    if (lr_disp[i] < 0)
      lr_disp[i] = nan;
    else
      lr_disp[i] *= -1.0;
  }
  
  char filename[] = "out_disp.tif";
  std::cout << "Writing " << filename << std::endl;
  iio_save_image_float((char*)filename, lr_disp, width, height);

#if 0
  char filename_lr[] = "lr_disp_pad.tif";
  std::cout << "Writing " << filename_lr << std::endl;
  iio_save_image_float((char*)filename_lr, lr_disp_pad, pad_width, height);

  char filename_rl[] = "rl_disp_pad.tif";
  std::cout << "Writing " << filename_rl << std::endl;
  iio_save_image_float((char*)filename_rl, rl_disp_pad, pad_width, height);
#endif
  
  // free memory
  free(l_img_pad);
  free(r_img_pad);
  free(lr_disp_pad);
  free(rl_disp_pad);
  free(lr_disp);
}

int main (int argc, char** argv) {

    process(argv[1], argv[2]);

  return 0;
}


