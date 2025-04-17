/*
Copyright 2011. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

This file is part of libelas.
Authors: Julius Ziegler, Andreas Geiger

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

#include <stdio.h>
#include <string.h>
#include <cassert>

#include "filter.h"

// define fixed-width datatypes for Visual Studio projects
#ifndef _MSC_VER
  #include <stdint.h>
#else
  typedef __int8            int8_t;
  typedef __int16           int16_t;
  typedef __int32           int32_t;
  typedef __int64           int64_t;
  typedef unsigned __int8   uint8_t;
  typedef unsigned __int16  uint16_t;
  typedef unsigned __int32  uint32_t;
  typedef unsigned __int64  uint64_t;
#endif

// fast filters: implements 3x3 and 5x5 sobel filters and 
//               5x5 blob and corner filters based on SSE2/3 instructions
namespace filter {
  
  // private namespace, public user functions at the bottom of this file
  namespace detail {
    void integral_image( const uint8_t* in, int32_t* out, int w, int h ) {
      int32_t* out_top = out;
      const uint8_t* line_end = in + w;
      const uint8_t* in_end   = in + w*h;
      int32_t line_sum = 0;
      for( ; in != line_end; in++, out++ ) {
        line_sum += *in;
        *out = line_sum;
      }
      for( ; in != in_end; ) {
        int32_t line_sum = 0;
        const uint8_t* line_end = in + w;
        for( ; in != line_end; in++, out++, out_top++ ) {
          line_sum += *in;
          *out = *out_top + line_sum;
        }
      }
    }
    
    void unpack_8bit_to_16bit( const __m128i a, __m128i& b0, __m128i& b1 ) {
      __m128i zero = _mm_setzero_si128();
      b0 = _mm_unpacklo_epi8( a, zero );
      b1 = _mm_unpackhi_epi8( a, zero );
    }
    
    void pack_16bit_to_8bit_saturate( const __m128i a0, const __m128i a1, __m128i& b ) {
      b = _mm_packus_epi16( a0, a1 );
    }
    
    // convolve image with a (1,2,1) row vector. Result is accumulated into output.
    // This one works on 16bit input and 8bit output.
    // output is scaled by 1/4, then clamped to [-128,128], and finally shifted to [0,255].
    void convolve_121_row_3x3_16bit( const int16_t* in, uint8_t* out, int w, int h ) {
      assert( w % 16 == 0 && "width must be multiple of 16!" );
      const __m128i* i0 = (const __m128i*)(in);
      const int16_t* i1 = in+1;
      const int16_t* i2 = in+2;
      uint8_t* result   = out + 1;
      const int16_t* const end_input = in + w*h;
      const size_t blocked_loops = (w*h-2)/16;
      __m128i offs = _mm_set1_epi16( 128 );
      for( size_t i=0; i != blocked_loops; i++ ) {
        __m128i result_register_lo;
        __m128i result_register_hi;
        __m128i i1_register;
        __m128i i2_register;
        
        i1_register        = _mm_loadu_si128( (__m128i*)( i1 ) );
        i2_register        = _mm_loadu_si128( (__m128i*)( i2 ) );
        result_register_lo = *i0;
        i1_register        = _mm_add_epi16( i1_register, i1_register );
        result_register_lo = _mm_add_epi16( i1_register, result_register_lo );
        result_register_lo = _mm_add_epi16( i2_register, result_register_lo );
        result_register_lo = _mm_srai_epi16( result_register_lo, 2 );
        result_register_lo = _mm_add_epi16( result_register_lo, offs );

        i0++;
        i1+=8;
        i2+=8;

        i1_register        = _mm_loadu_si128( (__m128i*)( i1 ) );
        i2_register        = _mm_loadu_si128( (__m128i*)( i2 ) );
        result_register_hi = *i0;
        i1_register        = _mm_add_epi16( i1_register, i1_register );
        result_register_hi = _mm_add_epi16( i1_register, result_register_hi );
        result_register_hi = _mm_add_epi16( i2_register, result_register_hi );
        result_register_hi = _mm_srai_epi16( result_register_hi, 2 );
        result_register_hi = _mm_add_epi16( result_register_hi, offs );

        i0++;
        i1+=8;
        i2+=8;

        pack_16bit_to_8bit_saturate( result_register_lo, result_register_hi, result_register_lo );
        _mm_storeu_si128( ((__m128i*)( result )), result_register_lo );
      
        result += 16;
      }
    }
    
    // convolve image with a (1,0,-1) row vector. Result is accumulated into output.
    // This one works on 16bit input and 8bit output.
    // output is scaled by 1/4, then clamped to [-128,128], and finally shifted to [0,255].
    void convolve_101_row_3x3_16bit( const int16_t* in, uint8_t* out, int w, int h ) {
      assert( w % 16 == 0 && "width must be multiple of 16!" );
      const __m128i*  i0 = (const __m128i*)(in);
      const int16_t* 	i2 = in+2;
      uint8_t* result    = out + 1;
      const int16_t* const end_input = in + w*h;
      const size_t blocked_loops = (w*h-2)/16;
      __m128i offs = _mm_set1_epi16( 128 );
      for( size_t i=0; i != blocked_loops; i++ ) {
        __m128i result_register_lo;
        __m128i result_register_hi;
        __m128i i2_register;

        i2_register = _mm_loadu_si128( (__m128i*)( i2 ) );
        result_register_lo  = *i0;
        result_register_lo  = _mm_sub_epi16( result_register_lo, i2_register );
        result_register_lo  = _mm_srai_epi16( result_register_lo, 2 );
        result_register_lo  = _mm_add_epi16( result_register_lo, offs );
 
        i0 += 1;
        i2 += 8;
        
        i2_register = _mm_loadu_si128( (__m128i*)( i2 ) );
        result_register_hi  = *i0;
        result_register_hi  = _mm_sub_epi16( result_register_hi, i2_register );
        result_register_hi  = _mm_srai_epi16( result_register_hi, 2 );
        result_register_hi  = _mm_add_epi16( result_register_hi, offs );

        i0 += 1;
        i2 += 8;
        
        pack_16bit_to_8bit_saturate( result_register_lo, result_register_hi, result_register_lo );
        _mm_storeu_si128( ((__m128i*)( result )), result_register_lo );

        result += 16;
      }

      for( ; i2 < end_input; i2++, result++) {
        *result = ((*(i2-2) - *i2)>>2)+128;
      }
    }
    
    void convolve_cols_3x3( const unsigned char* in, int16_t* out_v, int16_t* out_h, int w, int h ) {
      using namespace std;
      assert( w % 16 == 0 && "width must be multiple of 16!" );
      const int w_chunk  = w/16;
      __m128i* 	i0       = (__m128i*)( in );
      __m128i* 	i1       = (__m128i*)( in ) + w_chunk*1;
      __m128i* 	i2       = (__m128i*)( in ) + w_chunk*2;
      __m128i* result_h  = (__m128i*)( out_h ) + 2*w_chunk;
      __m128i* result_v  = (__m128i*)( out_v ) + 2*w_chunk;
      __m128i* end_input = (__m128i*)( in ) + w_chunk*h;
      for( ; i2 != end_input; i0++, i1++, i2++, result_v+=2, result_h+=2 ) {
        *result_h     = _mm_setzero_si128();
        *(result_h+1) = _mm_setzero_si128();
        *result_v     = _mm_setzero_si128();
        *(result_v+1) = _mm_setzero_si128();
        __m128i ilo, ihi;
        unpack_8bit_to_16bit( *i0, ihi, ilo ); 
        unpack_8bit_to_16bit( *i0, ihi, ilo );
        *result_h     = _mm_add_epi16( ihi, *result_h );
        *(result_h+1) = _mm_add_epi16( ilo, *(result_h+1) );
        *result_v     = _mm_add_epi16( *result_v, ihi );
        *(result_v+1) = _mm_add_epi16( *(result_v+1), ilo );
        unpack_8bit_to_16bit( *i1, ihi, ilo );
        *result_v     = _mm_add_epi16( *result_v, ihi );
        *(result_v+1) = _mm_add_epi16( *(result_v+1), ilo );
        *result_v     = _mm_add_epi16( *result_v, ihi );
        *(result_v+1) = _mm_add_epi16( *(result_v+1), ilo );
        unpack_8bit_to_16bit( *i2, ihi, ilo );
        *result_h     = _mm_sub_epi16( *result_h, ihi );
        *(result_h+1) = _mm_sub_epi16( *(result_h+1), ilo );
        *result_v     = _mm_add_epi16( *result_v, ihi );
        *(result_v+1) = _mm_add_epi16( *(result_v+1), ilo );
      }
    }
  };
  
  void sobel3x3( const uint8_t* in, uint8_t* out_v, uint8_t* out_h, int w, int h ) {
    int16_t* temp_h = (int16_t*)( _mm_malloc( w*h*sizeof( int16_t ), 16 ) );
    int16_t* temp_v = (int16_t*)( _mm_malloc( w*h*sizeof( int16_t ), 16 ) );    
    detail::convolve_cols_3x3( in, temp_v, temp_h, w, h );
    detail::convolve_101_row_3x3_16bit( temp_v, out_v, w, h );
    detail::convolve_121_row_3x3_16bit( temp_h, out_h, w, h );
    _mm_free( temp_h );
    _mm_free( temp_v );
  }
  
};
