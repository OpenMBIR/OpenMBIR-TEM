/* ============================================================================
 * Copyright (c) 2011, Michael A. Jackson (BlueQuartz Software)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of Michael A. Jackson nor the names of its contributors may
 * be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef MRCHEADER_H_
#define MRCHEADER_H_


#define FREE_FEI_HEADERS(feiHeaders)\
if (NULL != feiHeaders) {\
  free(feiHeaders);\
  feiHeaders = NULL;\
}\


/**
 * @brief this was take from http://www.biochem.mpg.de/doc_tom/index.html using the
 * tom_mrcfeistack2emseries code.
 */
typedef struct
{
    float a_tilt;  // float: 4 bytes
    float b_tilt;  // float: 4 bytes
    float x_stage; // float: 4 bytes
    float y_stage;  // float: 4 bytes
    float z_stage;  // float: 4 bytes
    float x_shift;  // float: 4 bytes
    float y_shift;  // float: 4 bytes
    float defocus;  // float: 4 bytes
    float exp_time; // float: 4 bytes
    float mean_int; // float: 4 bytes
    float tiltaxis; // float: 4 bytes
    float pixelsize;  // float: 4 bytes
    float magnification;    // float: 4 bytes
    float voltage; // float 4 bytes
    unsigned char unused[72]; // 72 bytes not used.
} FEIHeader;



/**
 * @brief This spec was taken from http://bio3d.colorado.edu/imod/doc/mrc_format.txt
 * and we are going to assume an IMOD version of 2.6.20 and above:
 */
typedef struct
{
/*    The MRC header. length 1024 bytes

    OFFSET  SIZE DATA    NAME Description */

/*      0 000  4 */ int     nx;  //Number of Columns
/*      4 004  4 */ int     ny;   //     Number of Rows
/*      8 010  4 */ int     nz;   //     Number of Sections.

/*      12 014  4 */ int     mode;     /* Types of pixel in image.  Values used by IMOD:
                       0 = unsigned or signed bytes depending on flag in imodStamp
                                 only unsigned bytes before IMOD 4.2.23
                       1 = signed short integers (16 bits)
                       2 = float
                       3 = short * 2, (used for complex data)
                       4 = float * 2, (used for complex data)
                       6 = unsigned 16-bit integers (non-standard)
                      16 = unsigned char * 3 (for rgb data, non-standard)
*/
/*     16 020  4 */ int     nxstart; //    Starting point of sub image (not used in IMOD)
/*     20 024  4 */ int     nystart;//
/*     24 030  4 */ int     nzstart;//

/*     28 034  4 */ int     mx;        // Grid size in X, Y, and Z
/*     32 040  4 */ int     my;//
/*     36 044  4 */ int     mz;//

/*     40 050  4 */ float   xlen;   //    Cell size; pixel spacing = xlen/mx, ylen/my, zlen/mz
/*     44 054  4 */ float   ylen;//
/*     48 060  4 */ float   zlen;//

/*     52 064  4 */ float   alpha;  //    cell angles - ignored by IMOD
/*     56 070  4 */ float   beta;//
/*     60 074  4 */ float   gamma;//

                  /*     These need to be set to 1, 2, and 3 for pixel spacing
                                       to be interpreted correctly */
/*     64 100  4 */ int     mapc; //      map column  1=x,2=y,3=z.
/*     68 104  4 */ int     mapr; //      map row     1=x,2=y,3=z.
/*     72 110  4 */ int     maps;  //     map section 1=x,2=y,3=z.

                                   //  These need to be set for proper scaling of data
/*     76 114  4 */ float   amin; //      Minimum pixel value.
/*     80 120  4 */ float   amax; //      Maximum pixel value.
/*     84 124  4 */ float   amean;  //    Mean pixel value.

/*     88 130  2 */ short   ispg;   //    image type
/*     90 132  2 */ short   nsymbt; //    space group number
/*     92 134  4 */ int     next;    //   number of bytes in extended header
/*     96 140  2 */ short   creatid; //   used to be an ID number, is 0 as of IMOD 4.2.23
/*     98 142  30 */ char  extra_data[30]; // (not used, first two bytes should be 0)

                          /*           These two values specify the structure of data in the
                                     extended header; their meaning depend on whether the
                                     extended header has the Agard format, a series of
                                     4-byte integers then real numbers, or has data
                                     produced by SerialEM, a series of short integers.
                                     SerialEM stores a float as two shorts, s1 and s2, by:
                                       value = (sign of s1)*(|s1|*256 + (|s2| modulo 256))
                                          * 2**((sign of s2) * (|s2|/256))   */
/*    128 200  2 */ short   nint;   //    Number of integers per section (Agard format) or
                                    // number of bytes per section (SerialEM format)
/*    130 202  2 */ short   nreal;  /*    Number of reals per section (Agard format) or bit
                                     flags for which types of short data (SerialEM format):
                                     1 = tilt angle * 100  (2 bytes)
                                     2 = piece coordinates for montage  (6 bytes)
                                     4 = Stage position * 25    (4 bytes)
                                     8 = Magnification / 100 (2 bytes)
                                     16 = Intensity * 25000  (2 bytes)
                                     32 = Exposure dose in e-/A2, a float in 4 bytes
                                     128, 512: Reserved for 4-byte items
                                     64, 256, 1024: Reserved for 2-byte items
                                     If the number of bytes implied by these flags does
                                     not add up to the value in nint, then nint and nreal
                                     are interpreted as ints and reals per section
*/
/*    132 204  20 */ char     extra_data_2[20]; // (not used)

/*    152 230  4 */ int     imodStamp; //   1146047817 indicates that file was created by IMOD or
                                     //  other software that uses bit flags in the following field
/*    156 234  4 */ int     imodFlags; //   Bit flags:
                                       // 1 = bytes are stored as signed

                              //    Explanation of type of data.
/*    160 240  2 */ short   idtype; // ( 0 = mono, 1 = tilt, 2 = tilts, 3 = lina, 4 = lins)
/*    162 242  2 */ short   lens;//
/*    164 244  2 */ short   nd1; //for idtype = 1, nd1 = axis (1, 2, or 3)
/*    166 246  2 */ short   nd2;//
/*    168 250  2 */ short   vd1;    //                   vd1 = 100. * tilt increment
/*    170 252  2 */ short   vd2;   //                    vd2 = 100. * starting angle

     /*                   Current angles are used to rotate a model to match a
                                    new rotated image.  The three values in each set are
                                    rotations about X, Y, and Z axes, applied in the order
                                    Z, Y, X. */
/*    172 254  24 */ float   tiltangles[6]; // 0,1,2 = original:  3,4,5 = current
/*
                                   The image origin is the location of the origin of the
                                   coordinate system relative to the first pixel in the
                                   file.  It is in pixel spacing units rather than in
                                   pixels.  If an original volume has an origin of 0, a
                                   subvolume should have negative origin values.
                                   */

/* NEW-STYLE MRC image2000 HEADER - IMOD 2.6.20 and above: */
/* 196 304  4 */  float   xorg;   //   Origin of image
/* 200 310  4 */ float   yorg;//
/* 204 314  4 */ float   zorg;//

/* 208 320  4 */ char    cmap[4];     // Contains "MAP "
/* 212 324  4 */ char    stamp[4];   //  First two bytes have 17 and 17 for big-endian or
                                     //   68 and 65 for little-endian
/* 216 330  4 */ float   rms;      // RMS deviation of densities from mean density

//ALL HEADERS:
/* 220 334  4 */ int     nLabels;  // Number of labels with useful data.
/* 224 340  800 */ char   labels[10][80]; //  10 labels of 80 charactors.

/*
------------------------------------------------------------------------

Offsets are given as bytes in decimal and octal.
Total size of header is 1024 bytes plus the size of the extended header.
*/

/**
 * @brief Looks like on FEI instruments there are values that can be parsed
 * from the extended header. The number of entries depends on the number of tilts
 * that were saved to the file (nz).
 */

  FEIHeader*   feiHeaders;

/*
Image data follows with the origin in the lower left corner,
looking down on the volume.

The size of the image is nx * ny * nz * (mode data size).
*/

} MRCHeader;



#endif /* MRCHEADER_H_ */
