This is a fork of ``libelas``, from

    https://github.com/kou1okada/libelas

It is released under the same terms as the original.

The changes are at the interface level, the resulting tool is able to
parse input options and read and write TIF files with float values. It
can handle negative disparities and command line options. It changes
the sign of the disparity before writing it back to disk, applies some
outlier filtering to the disparity, and sets invalid disparities to
NaN.

The program assumes that the input images have float pixels with values 
between 0 and 1. Those will be scaled to 0 to 255 and cast to uint8 before 
invoking libelas. Use the option ``-scale val`` to change the value of the 
multiplier.

The facility for handling TIF files is borrowed from
https://github.com/gfacciol/mgm.

Usage: 

    elas <options> left_image.tif right_image.tif output_disparity.tif

Invoking this program with the correct inputs and the opiton `-verbose
1`, will print the values of the parameters being used by this
program. Those can be overridden from the command line, for example
as:

    elas -ipol_gap_width 100 -postprocess_only_left 0 <other options>

Command line options
====================

-disp_min (default = 0):
    Minimum disparity (ASP estimates this unless the user overrides it).

-disp_max (default = 255):
    Maximum disparity (ASP estimates this unless the user overrides it).

-support_threshold (default = 0.85):
    Maximum uniqueness ratio (best vs. second-best support match).

-support_texture (default = 10):
    Minimum texture for support points.

-candidate_stepsize (default = 5):
    Step size of regular grid on which support points are matched.

-incon_window_size (default = 5):
    Window size of inconsistent support point check.

-incon_threshold (default = 5):
    Disparity similarity threshold for support point to be considered
    consistent.

-incon_min_support (default = 5):
    Minimum number of consistent support points.

-add_corners (default = 0):
    Add support points at image corners with nearest neighbor
    disparities.

-grid_size (default = 20):
    Size of neighborhood for additional support point extrapolation.

-beta (default = 0.02):
    Image likelihood parameter.

-gamma (default = 3):
    Prior constant.

-sigma (default = 1):
    Prior sigma.

-sradius (default = 2):
    Prior sigma radius.

-match_texture (default = 1):
    Minimum texture for dense matching.

-lr_threshold (default = 2):
    Disparity threshold for left-right consistency check.

-speckle_sim_threshold (default = 1):
    Similarity threshold for speckle segmentation.

-speckle_size (default = 200):
    Speckles larger than this get removed.

-ipol_gap_width (default = 3):
    Fill holes in disparity of height and width at most this value.

-filter_median (default = 0):
    If non-zero, use an approximate median filter.

-filter_adaptive_mean (default = 1):
    If non-zero, use an approximate adaptive mean filter.

-postprocess_only_left (default = 0):
    If non-zero, saves time by not postprocessing the right image.

-verbose (default = 0):
    If non-zero, print some information about the values of the
    options being used, as well as what the input and output files
    are.

-scale (default = 255.0):

   Multiply the input pixels by this value before casting to uint8. Hence the 
   default assumes that the pixels are between 0 and 1.
   
Original documentation is below.
===============================

####################################################################################
# Copyright 2011. All rights reserved.                                             #
# Institute of Measurement and Control Systems                                     #
# Karlsruhe Institute of Technology, Germany                                       #
#                                                                                  #
# This file is part of libelas.                                                    #
# Authors:  Andreas Geiger                                                         #
#           Please send any bugreports to geiger@kit.edu                           #
#                                                                                  #
# libelas is free software; you can redistribute it and/or modify it under the     #
# terms of the GNU General Public License as published by the Free Software        #
# Foundation; either version 3 of the License, or any later version.               #
#                                                                                  #
# libelas is distributed in the hope that it will be useful, but WITHOUT ANY       #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.         #
#                                                                                  #
# You should have received a copy of the GNU General Public License along with     #
# libelas; if not, write to the Free Software Foundation, Inc., 51 Franklin        #
# Street, Fifth Floor, Boston, MA 02110-1301, USA                                  #
####################################################################################

+++++++++++++++++++++++++++++++++++
+          INTRODUCTION           +
+++++++++++++++++++++++++++++++++++

Libelas (LIBrary for Efficient LArge-scale Stereo matching) is a cross-platfrom C++
library with MATLAB wrappers for computing disparity maps of large images. Input
is a rectified greyscale stereo image pair of same size. Output are the corresponding
disparity maps.

If you distribute a software that uses libelas, you have to distribute it under GPL
with the source code. Another option is to contact us to purchase a commercial license.

If you find this software useful or if you use this software for your research,
we would be happy if you cite the following related publication:

@INPROCEEDINGS{Geiger10,
 author = {Andreas Geiger and Martin Roser and Raquel Urtasun},
 title = {Efficient Large-Scale Stereo Matching},
 booktitle = {Asian Conference on Computer Vision},
 year = {2010},
 month = {November},
 address = {Queenstown, New Zealand}
}

+++++++++++++++++++++++++++++++++++
+    COMPILING MATLAB WRAPPERS    +
+++++++++++++++++++++++++++++++++++

If you want to use libelas directly from MATLAB you can easily do this by using
the MATLAB wrappers provided. They also include some demo files for testing your
configuration. First, configure your MATLAB MEX C++ compiler, if it is not yet
configured (mex -setup). Under Linux you might use g++, under Windows I compiled
it successfully with the Microsoft Visual Studio Express 2008 compilers.

1) Change to the libelas/matlab directory
2) After running 'make.m' you should have a MEX file called 'elasMex'
3) Now try to run 'demo.m' which opens a GUI and shows you
   some results on the included test images

+++++++++++++++++++++++++++++++++++
+     BUILDING A C++ LIBRARY      +
+++++++++++++++++++++++++++++++++++

Prerequisites needed for compiling libelas using c++:
- CMake (available at: http://www.cmake.org/)

Linux:

1) Move to libelas root directory
2) Type 'cmake .'
3) Type 'make'
4) Run './elas demo' => computes disparity maps for images from the 'img' directory

Windows:

1) Start CMake GUI
2) Set directories to elas root directory
3) Run configure, configure and generate
4) Open the resulting Visual Studio solution with Visual Studio
5) Switch to 'Release' mode and build all
6) Move 'elas.exe' from libelas/Release to libelas
7) Open a console and navigate to libelas root directory
8) Run 'elas.exe demo' => computes disparity maps in the img directory

For more information on CMake, have a look at the CMake documentation.

For more information on the usage of the library, have a look into the MATLAB wrappers and
into the documentation of the header elas.h.

Please send any feedback and bugreports to geiger@kit.edu
Andreas Geiger

