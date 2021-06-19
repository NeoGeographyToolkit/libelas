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
#include <map>
#include <tiffio.h>
#include <cmath>
#include "elas.h"
#include "image.h"

extern "C" {
#include "iio.h"
}

using namespace std;

// These ugly macros saves a lot of typing

#define SET_VAL(val)               \
  it = options.find("-"#val);      \
  if (it != options.end())         \
    params.val = it->second;

#define PRINT_VAL(val)             \
  std::cout << #val << " = " << params.val << std::endl;

// Fill the parameters from any options the user specified.  A lot of
// repetitive code in this function.
void fill_params_from_options(std::map<std::string, float> const & options,
                              Elas::parameters & params) {

  std::map<std::string, float>::const_iterator it;
  SET_VAL(disp_min);
  SET_VAL(disp_max);
  SET_VAL(support_threshold);
  SET_VAL(support_texture);
  SET_VAL(candidate_stepsize);
  SET_VAL(incon_window_size);
  SET_VAL(incon_threshold);
  SET_VAL(incon_min_support);
  SET_VAL(add_corners);
  SET_VAL(grid_size);
  SET_VAL(beta);
  SET_VAL(gamma);
  SET_VAL(sigma);
  SET_VAL(sradius);
  SET_VAL(match_texture);
  SET_VAL(lr_threshold);
  SET_VAL(speckle_sim_threshold);
  SET_VAL(speckle_size);
  SET_VAL(ipol_gap_width);
  SET_VAL(filter_median);
  SET_VAL(filter_adaptive_mean);
  SET_VAL(postprocess_only_left);
}

void print_param_vals(Elas::parameters const& params) {
  std::cout << "Invoking LIBELAS with the following values:\n";
  PRINT_VAL(disp_min);
  PRINT_VAL(disp_max);
  PRINT_VAL(support_threshold);
  PRINT_VAL(support_texture);
  PRINT_VAL(candidate_stepsize);
  PRINT_VAL(incon_window_size);
  PRINT_VAL(incon_threshold);
  PRINT_VAL(incon_min_support);
  PRINT_VAL(add_corners);
  PRINT_VAL(grid_size);
  PRINT_VAL(beta);
  PRINT_VAL(gamma);
  PRINT_VAL(sigma);
  PRINT_VAL(sradius);
  PRINT_VAL(match_texture);
  PRINT_VAL(lr_threshold);
  PRINT_VAL(speckle_sim_threshold);
  PRINT_VAL(speckle_size);
  PRINT_VAL(ipol_gap_width);
  PRINT_VAL(filter_median);
  PRINT_VAL(filter_adaptive_mean);
  PRINT_VAL(postprocess_only_left);
}

#undef SET_VAL
#undef PRINT_VAL

int process(std::map<std::string, float> const& options,
            std::string const& left_image, std::string const& right_image,
            std::string const& out_disp) {
  
  std::cout << "Processing: " << left_image << ", " << right_image << std::endl;

  // TODO(oalexan1): Need to think more here.
  Elas::parameters params;
  fill_params_from_options(options, params);

  std::cout << "Input disp_min and disp_max: "
            << params.disp_min << ' ' << params.disp_max << std::endl;
  
  // Note how we adjust disp_min and disp_max, as libelas likes only a positive disparity
  int pad = std::max(-params.disp_min, 0);
  params.disp_min = 0;
  params.disp_max = params.disp_max + pad;
  std::cout << "Adjust for LIBELAS only being able to handle positive disparities." << std::endl;
  std::cout << "Pad the left image on the left with " << pad << " columns of zeros."<< std::endl;

  int width, height; 
  float * limg = iio_read_image_float(left_image.c_str(), &width, &height);

  int rwidth, rheight; 
  float * rimg = iio_read_image_float(right_image.c_str(), &rwidth, &rheight);

  // check for correct size
  if (width <= 0 || height <= 0 || rwidth <= 0 || rheight <= 0 ||
      width != rwidth || height != rheight) {
    std::cout << "ERROR: Images must be of same size, but got: " << std::endl;
    std::cout << "  left_image:  " << width <<  " x " << height  << std::endl;
    std::cout << "  right_image: " << rwidth <<  " x " << rheight << std::endl;
    return 1; 
  }

  int pad_width = width + pad;

  // TODO(oalexan1): Make these into functions.
  
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
  Elas elas(params);
  print_param_vals(params);
  elas.process(l_img_pad, r_img_pad, lr_disp_pad, rl_disp_pad, dims);

  // When the disparities have a big jump, the ones after the jump are outliers.
  // TODO(oalexan1): This is fragile.
  // TODO(oalexan1): Make this into a function.
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
  
  // Remove the padding by cropping the disparity
  count = 0;
  float* lr_disp = (float*)malloc(width * height * sizeof(float));
  for (int ih = 0; ih < height; ih++) {
    for (int iw = 0; iw < width; iw++) {
      lr_disp[count] = lr_disp_pad[ih * pad_width + iw + pad];
      count++;
    }
  }
  
  // Negative disparities are set to NaN.  Subtract the padding. Flip
  // the sign of the disparity to follow the ASP convention. Ignore
  // the rl disparity.
  float nan = std::numeric_limits<float>::quiet_NaN();
  for (int32_t i = 0; i < width * height; i++) {
    if (lr_disp[i] < 0)
      lr_disp[i] = nan;
    else
      lr_disp[i] = -(lr_disp[i] - pad);
  }
  
  std::cout << "Writing " << out_disp << std::endl;
  iio_save_image_float((char*)out_disp.c_str(), lr_disp, width, height);

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

  return 0;
}

// Parse the input arguments. Assume that each option starts with a dash,
// and its value is a float. The rest of the arguments are files.
int parse_args(int argc, char** argv,
                std::map<std::string, float> & options,
                std::vector<std::string> & files) {

  // Wipe the outputs first
  options.clear();
  files.clear();

  for (int it = 1; it < argc; it++) {

    if (strlen(argv[it]) > 0 && argv[it][0] == '-') {
      // This is an option
      if (it == argc - 1) {
        std::cout << "ERROR: Found an option without a value." << std::endl;
        return 1;
      }
      options[std::string(argv[it])] = atof(argv[it+1]);

      // Jump over the value we just parsed
      it++;
    } else {
      files.push_back(argv[it]);
    }
    
  }
  return 0;
}
  
int main (int argc, char** argv) {
  
  // Parse the options 
  std::map<std::string, float> options;
  std::vector<std::string> files;
  if (parse_args(argc, argv, options, files) != 0)
    return 1;
  
  if (files.size() < 3) {
    std::cout << "Usage: elas <options> left_image.tif right_image.tif output_disp.tif\n";
    return 1;
  }
  
  std::string left_image, right_image, out_disp;
  left_image  = files[0];
  right_image = files[1];
  out_disp    = files[2];

  std::cout << "Reading " << left_image << ' ' << right_image << ' ' << out_disp << std::endl;

  return process(options, left_image, right_image, out_disp);
}


