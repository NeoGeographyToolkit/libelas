#ifndef _IIO_H
#define _IIO_H

#include <stdint.h>
#include <stdbool.h>


// Original is from https://github.com/gfacciol/mgm/blob/master/iio/iio.h

//
//
//
// IIO: "Image Input-Output"
//
// A library for reading and writing small images
//
//
//


///////////////////////////////
// HIGH-LEVEL API FUNCTIONS  //
// (using only C constructs) //
///////////////////////////////



//
// basic float API for 2D images (returns a freeable pointer)
//

float *iio_read_image_float(const char *fname, int *w, int *h);
// x[i+j*w]

float *iio_read_image_float_vec(const char *fname, int *w, int *h, int *pd);
// x[(i + j*w)*pd + l]

float *iio_read_image_float_rgb(const char *fname, int *w, int *h);

float *iio_read_image_float_split(const char *fname, int *w, int *h, int *pd);
// x[w*h*l + i + j*w]

//
// convenience float API for 2D images (also returns a freeable pointer)
//

float **iio_read_image_float_matrix(const char *fname, int *w, int *h);
// x[j][i]

float (**iio_read_image_float_matrix_2vec(const char *fnam, int *w, int *h))[2];
float (**iio_read_image_float_matrix_3vec(const char *fnam, int *w, int *h))[3];
float (**iio_read_image_float_matrix_4vec(const char *fnam, int *w, int *h))[4];
float (**iio_read_image_float_matrix_rgb(const char *fnam, int *w, int *h))[3];
float (**iio_read_image_float_matrix_rgba(const char *fnam, int *w, int *h))[4];
// x[j][i][channel]
// (The "rgb" and "rgba" functions may re-order the channels according to file
// metadata.  The "vec" functions produce the data in the same order as is
// stored.)
void *iio_read_image_float_matrix_vec(const char *fnam, int *w, int *h, int *pd);



//
// So far, we have seen the functions for reading multi-channel "2D" images,
// whose samples are of type "FLOAT", from a "NAMED FILE".  The rest of the
// API consists in variations of these functions, changing the quoted options
// of the previous sentence in all possible ways.  For instance, there are
// the following variations of the function "iio_read_image_float":
//
// 	// read 3D image into floats, from a named file
// 	float *iio_read_3d_image_float(char *fname, int *w, int *h, int *d);
//
//	// read 2D image into bytes, from a named file
// 	uint8_t *iio_read_image_uint8(char *fname, int *w, int *h);
//
// 	// read 2D image into floats, from an open stream
// 	float *iio_read_image_float_f(FILE *f, int *w, int *h);
//
// 	// read 2D image into floats, from memory
// 	float *iio_read_image_float_m(void *f, size_t s, int *w, int *h);
//
//	// general forms of these functions
//	TYPE *iio_read{,_3d,_4d,_nd}_image_TYPE{,_f,_m}(...);
//	TYPE **iio_read{,_3d,_4d,_nd}_image_TYPE_matrix{,_f,_m}(...);
//	TYPE (**iio_read{,_3d,_4d,_nd}_image_TYPE_matrix_2vec{,_f,_m}(...))[2];
//	TYPE (**iio_read{,_3d,_4d,_nd}_image_TYPE_matrix_3vec{,_f,_m}(...))[3];
//	TYPE (**iio_read{,_3d,_4d,_nd}_image_TYPE_matrix_4vec{,_f,_m}(...))[4];
//	TYPE (**iio_read{,_3d,_4d,_nd}_image_TYPE_matrix_rgb{,_f,_m}(...))[3];
//	TYPE (**iio_read{,_3d,_4d,_nd}_image_TYPE_matrix_rgba{,_f,_m}(...))[4];
//
//
double *iio_read_image_double(const char *fname, int *w, int *h);
double *iio_read_image_double_vec(const char *fname, int *w, int *h, int *pd);

int *iio_read_image_int(const char *fname, int *w, int *h);

// All these functions are boring  variations, and they are defined at the
// end of this file.  More interesting are the two following general
// functions:

void *iio_read_image_numbers_as_they_are_stored(char *fname, int *w, int *h,
		int *samples_per_pixel, int *sample_size,
		bool *ieeefp_samples, bool *signed_samples);

void *iio_read_image_numbers_in_requested_format(char *fname, int *w, int *h,
		int *samples_per_pixel, int requested_sample_size,
		bool requested_ieeefp, bool requested_sign);

// These two general functions have the usual versions for nD images and
// streams.  There exist also the following truly general functions, that
// read images of any dimension:

void *iio_read_nd_image_as_stored(char *fname,
		int *dimension, int *sizes,
		int *samples_per_pixel, int *sample_size,
		bool *ieeefp_samples, bool *signed_samples);

void *iio_read_nd_image_as_desired(char *fname,
		int *dimension, int *sizes,
		int *samples_per_pixel, int desired_sample_size,
		bool desired_ieeefp_samples, bool desired_signed_samples);






#ifdef UINT8_MAX

// basic byte API (returns a freeable pointer)
//

uint8_t *iio_read_image_uint8(const char *fname, int *w, int *h);
// x[i+j*w]

uint8_t *iio_read_image_uint8_vec(const char *fname, int *w, int *h, int *nc);
// x[(i + j*w)*nc + l]
//
uint8_t (*iio_read_image_uint8_rgb(const char *fnam, int *w, int *h))[3];


//
// convenience float API (also returns a freeable pointer)
//

uint8_t **iio_read_image_uint8_matrix(const char *fname, int *w, int *h);
// x[j][i]

uint8_t (**iio_read_image_uint8_matrix_2vec(const char *fnam, int *w, int *h))[2];
uint8_t (**iio_read_image_uint8_matrix_3vec(const char *fnam, int *w, int *h))[3];
uint8_t (**iio_read_image_uint8_matrix_4vec(const char *fnam, int *w, int *h))[4];
uint8_t (**iio_read_image_uint8_matrix_rgb(const char *fnam, int *w, int *h))[3];
uint8_t (**iio_read_image_uint8_matrix_rgba(const char *fnam, int *w, int *h))[4];
// x[j][i][channel]
// (The "rgb" and "rgba" functions may re-order the channels according to file
// metadata.  The "vec" functions produce the data in the same order as is
// stored.)
uint8_t ***iio_read_image_uint8_matrix_vec(const char *fnam, int *w, int *h, int *pd);

#endif//UINT8_MAX

#ifdef UINT16_MAX
uint16_t *iio_read_image_uint16_vec(const char *fname, int *w, int *h, int *pd);
#endif//UINT16_MAX

void iio_save_image_float_vec(char *filename, float *x, int w, int h, int pd);
void iio_save_image_float_split(char *filename, float *x, int w, int h, int pd);
void iio_save_image_double_vec(char *filename, double *x, int w, int h, int pd);
void iio_save_image_float(char *filename, float *x, int w, int h);
void iio_save_image_double(char *filename, double *x, int w, int h);
void iio_save_image_int(char *filename, int *x, int w, int h);
void iio_save_image_int_vec(char *filename, int *x, int w, int h, int pd);
void iio_save_image_uint8_vec(char *filename, uint8_t *x, int w, int h, int pd);
void iio_save_image_uint16_vec(char *filename, uint16_t *x, int w, int h, int pd);
void iio_save_image_uint8_matrix_rgb(char *f, unsigned char (**x)[3], int w, int h);
void iio_save_image_uint8_matrix(char *f, unsigned char **x, int w, int h);

// SAVING FORMATS:
// (w, h; 1 uint8) => pgm, png, tiff, pfm
// (w, h; 3 uint8) => ppm, png, tiff, pfm
// (w, h; 1 uint16) => pgm, png, tiff, pfm
// (w, h; 3 uint16) => ppm, png, tiff, pfm
// (w, h; 1 float) => tiff, pfm
// (w, h; 3 float) => tiff, pfm

#define IIO_MAX_DIMENSION 20

#endif//_IIO_H
