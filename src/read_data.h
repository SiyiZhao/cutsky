/*******************************************************************************
* read_data.h: this file is part of the cutsky program.

* cutsky: cutsky catalogue generator.

* Github repository:
        https://github.com/cheng-zhao/cutsky

* Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

*******************************************************************************/

#ifndef __READ_DATA_H__
#define __READ_DATA_H__

#include <stdio.h>

/*============================================================================*\
                        Functions for reading data files
\*============================================================================*/

/******************************************************************************
Function `read_ascii_twocol`:
  Read the first two columns of an ASCII file as double arrays.
Arguments:
  * `fname`:    filename of the input catalog;
  * `x`:        array for the first column;
  * `y`:        array for the second column;
  * `num`:      number of lines read successfully.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_ascii_twocol(const char *fname, double **x, double **y, size_t *num);


#ifdef WITH_CFITSIO

/******************************************************************************
Function `read_filelist`:
  Read filenames from a list.
Arguments:
  * `fname`:    filename of the input list;
  * `list`:     list of filenames read from file;
  * `num`:      number of lines read successfully.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_filelist(const char *fname, char ***list, int *num);

#endif

#endif
