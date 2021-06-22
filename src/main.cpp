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

// Parse the input arguments. Assume that each option starts with a
// dash, and it has only one value, which is a float. The rest of the
// arguments are files.
int parse_args(int argc, char** argv,
               std::map<std::string, float> & options,
               std::vector<std::string> & files) {
  
  // Wipe the outputs first
  options.clear();
  files.clear();

  for (int it = 1; it < argc; it++) {

    // Sanity check. There must be no double-dashes.
    if (strlen(argv[it]) > 1 && argv[it][0] == '-' && argv[it][1] == '-') {
      std::cout << "ERROR: Only single-dash options are allowed. Got: " << argv[it] << std::endl;
      return 1;
    }
    
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

// Copy an image as a block in a wider image with the same height,
// scale the pixels, and transform to uint8. Assume that the images
// have been allocated by now.
void insert_scale_block(float const* img, uint8_t* img_pad,
                        double scale, int width, int height,
                        int padded_width, int pad) {

  for (int i = 0; i < padded_width * height; i++) 
    img_pad[i] = 0;
  
  int count = 0;
  for (int ih = 0; ih < height; ih++) {
    for (int iw = 0; iw < width; iw++) {
      
      float val = img[count];
      if (std::isnan(val)) 
        val = 0.0;
      
      val *= scale;
      val = round(val);
      if (val < 0) 
        val = 0;
      if (val > 255.0)
        val = 255.0;
      
      // Insert padding on the left 
      img_pad[ih * padded_width + iw + pad] = (uint8_t)val;
      
      count++;
    }
  }

}

// When the disparities have a big jump, the ones after the jump are
// outliers.  Remove those.
// TODO(oalexan1): This is fragile.
void filter_disparity(float * disp, int width, int height) {
  
  float* sorted_disp = (float*)malloc(width * height * sizeof(float));
  for (int i = 0; i < width * height; i++)
    sorted_disp[i] = disp[i];
  std::sort(sorted_disp, sorted_disp + width * height);

  float max_jump = -1.0, max_valid = -1.0;
  for (int i = 0; i < width * height - 1; i++) {
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

  // TODO(oalexan1): Need to think more here
  if (max_jump > 2) {
    for (int i = 0; i < width * height; i++)
      if (disp[i] > max_valid)
        disp[i] = -10.0; // invalidate the outlier disparity
  }
  
}
void savePGM(int width, int height, uint8_t* img, const char *name) {

  int uchar_max = 255;
  std::ofstream file(name, std::ios::out | std::ios::binary);
  file << "P5\n" << width << " " << height << "\n" << uchar_max << "\n";
  file.write((char *)img, width * height * sizeof(uint8_t));
}

int process(std::map<std::string, float> const& options,
            std::string const& left_image, std::string const& right_image,
            std::string const& out_disp) {
  
  Elas::parameters params;
  fill_params_from_options(options, params);

  bool verbose = false;
  auto it = options.find("-verbose");
  if (it != options.end() && it->second != 0) {
    verbose = true;
  }

  if (verbose) 
    std::cout << "Input disp_min and disp_max: "
              << params.disp_min << ' ' << params.disp_max << std::endl;

  // Multiply image pixels by this value
  double scale = 255.0;
  it = options.find("-scale");
  if (it != options.end() && it->second != 0) {
    scale = it->second;
  }

  if (verbose) 
    std::cout << "Multiplying input pixels by: " << scale << std::endl;
  
  // Note how we adjust disp_min and disp_max, as libelas likes only a positive disparity
  int pad = std::max(-params.disp_min, 0);
  params.disp_min = 0;
  params.disp_max = params.disp_max + pad;
  if (verbose) {
    std::cout << "Adjust for LIBELAS only being able to handle positive disparities." << std::endl;
    std::cout << "Pad the left image on the left with " << pad << " columns of zeros."<< std::endl;
    std::cout << "Reading: " << left_image << ' ' << right_image << std::endl;
  }
  
  int width, height; 
  float * l_img = iio_read_image_float(left_image.c_str(), &width, &height);

  int rwidth, rheight; 
  float * r_img = iio_read_image_float(right_image.c_str(), &rwidth, &rheight);

  // check for correct size
  if (width <= 0 || height <= 0 || rwidth <= 0 || rheight <= 0 ||
      width != rwidth || height != rheight) {
    std::cout << "ERROR: Images must be of same size, but got: " << std::endl;
    std::cout << "  left_image:  " << width <<  " x " << height  << std::endl;
    std::cout << "  right_image: " << rwidth <<  " x " << rheight << std::endl;
    return 1; 
  }

  int padded_width = width + pad;

  // Convert the image pixels from floats in [0, 1] to uint8_t in [0,
  // 255].  To ensure a positive disparity, pad the left image on the
  // left with enough zeros, and pad the right image on the right with
  // the same amount.

  // Process the left image
  uint8_t* l_img_pad = (uint8_t*)malloc(padded_width * height * sizeof(uint8_t));
  insert_scale_block(l_img, l_img_pad, scale, width, height, padded_width, pad);

  // Process the right image. Note how we pad right, so we insert no extra
  // columns on the left.
  uint8_t* r_img_pad = (uint8_t*)malloc(padded_width * height * sizeof(uint8_t));
  insert_scale_block(r_img, r_img_pad, scale, width, height, padded_width, 0);

#if 0
  // For debugging, save as tif
  int ch = 1;
  char l_img_pad_file[] = "l_img_pad.tif";
  std::cout << "Writing: " << l_img_pad_file << std::endl;
  iio_save_image_uint8_vec((char*)l_img_pad_file, l_img_pad, padded_width, height, ch);

  char r_img_pad_file[] = "r_img_pad.tif";
  std::cout << "Writing: " << r_img_pad_file << std::endl;
  iio_save_image_uint8_vec((char*)r_img_pad_file, r_img_pad, padded_width, height, ch);
#endif

#if 0
  // for debugging, save as pgm
  char l_img_pad_file[] = "l_img_pad.pgm";
  std::cout << "Writing: " << l_img_pad_file << std::endl;
  savePGM(padded_width, height, l_img_pad, l_img_pad_file);

  char r_img_pad_file[] = "r_img_pad.pgm";
  std::cout << "Writing: " << r_img_pad_file << std::endl;
  savePGM(padded_width, height, r_img_pad, r_img_pad_file);
#endif

  // allocate memory for disparity images
  float* lr_disp_pad = (float*)malloc(padded_width * height * sizeof(float));
  float* rl_disp_pad = (float*)malloc(padded_width * height * sizeof(float));
  
  // process
  const int32_t dims[3] = {padded_width, height, padded_width}; // bytes per line = width
  Elas elas(params);
  if (verbose) 
    print_param_vals(params);
  elas.process(l_img_pad, r_img_pad, lr_disp_pad, rl_disp_pad, dims);

  // When the disparities have a big jump, the ones after the jump are outliers.
  // TODO(oalexan1): This is fragile.
  filter_disparity(lr_disp_pad, padded_width, height);

  // Remove the padding by cropping the disparity
  int count = 0;
  float* lr_disp = (float*)malloc(width * height * sizeof(float));
  for (int ih = 0; ih < height; ih++) {
    for (int iw = 0; iw < width; iw++) {
      lr_disp[count] = lr_disp_pad[ih * padded_width + iw + pad];
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

  if (verbose) 
    std::cout << "Writing " << out_disp << std::endl;
  iio_save_image_float((char*)out_disp.c_str(), lr_disp, width, height);

#if 0
  char filename_lr[] = "lr_disp_pad.tif";
  std::cout << "Writing " << filename_lr << std::endl;
  iio_save_image_float((char*)filename_lr, lr_disp_pad, padded_width, height);

  char filename_rl[] = "rl_disp_pad.tif";
  std::cout << "Writing " << filename_rl << std::endl;
  iio_save_image_float((char*)filename_rl, rl_disp_pad, padded_width, height);
#endif
  
  // free memory
  free(l_img_pad);
  free(r_img_pad);
  free(lr_disp_pad);
  free(rl_disp_pad);
  free(lr_disp);

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

  return process(options, left_image, right_image, out_disp);
}


