/*******************************************************************************
* read_file.h: this file is part of the cutsky program.

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

#ifndef __READ_FILE_H__
#define __READ_FILE_H__

#include <stdio.h>

#ifdef WITH_CFITSIO
#include <fitsio.h>
#endif

/*============================================================================*\
                        Data structures for file reading
\*============================================================================*/

typedef struct {
  const char *fname;    /* name of the input file                        */
  FILE *fp;             /* pointer to the file stream                    */
  char *chunk;          /* buffer for reading file by chunks             */
  const char *line;     /* pointer to the current record                 */
  size_t *lines;        /* offsets of all lines in the chunk             */
  size_t used;          /* number of processed bytes in the buffer       */
  size_t rest;          /* number of bytes to be processed in the buffer */
  size_t size;          /* capacity of the buffer                        */
  size_t nline;         /* number of lines in the chunk                  */
  size_t maxline;       /* capacity of the array storing lines           */
} IFILE;

#ifdef WITH_CFITSIO
typedef struct {
  const char **fnames;  /* names of input files                          */
  int ninput;           /* number of input files                         */
  int current;          /* currently opened file                         */
  fitsfile *fp;         /* pointer to the FITS file stream               */
  int col[6];           /* column numbers of coordinates and velocities  */
  double *data[6];      /* coordinates and velocities read from file     */
  long ntotal;          /* number of objects in the current file         */
  long nread;           /* number of processed objects in the file       */
  long nstep;           /* optimal number of rows to be read at once     */
  long ndata;           /* number of reported objects                    */
  long nsave;           /* number of stored objects                      */
  long maxdata;         /* capacity of the arrays for reported data      */
} IFFILE;
#endif

/*============================================================================*\
                       Interfaces for ASCII file reading
\*============================================================================*/

/******************************************************************************
Function `input_init`:
  Initialise the interface for file reading.
Return:
  Address of the interface.
******************************************************************************/
IFILE *input_init(void);

/******************************************************************************
Function `input_destroy`:
  Deconstruct the interface for file reading.
Arguments:
  * `ifile`:    structure for file reading.
******************************************************************************/
void input_destroy(IFILE *ifile);

/******************************************************************************
Function `input_newfile`:
  Open a new file for reading.
Arguments:
  * `ifile`:    interface for file reading;
  * `fname`:    name of the file to be read from.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int input_newfile(IFILE *ifile, const char *fname);

/******************************************************************************
Function `input_readline`:
  Read a record (line) from the input file.
Arguments:
  * `ifile`:    interface for file reading.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int input_readline(IFILE *ifile);

/******************************************************************************
Function `input_readlines`:
  Read multiple records (lines) from the input file.
Arguments:
  * `ifile`:    interface for file reading;
  * `nline`:    number of lines to be read at once.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int input_readlines(IFILE *ifile, const size_t nline);


#ifdef WITH_CFITSIO

/*============================================================================*\
                        Interfaces for FITS file reading
\*============================================================================*/

/******************************************************************************
Function `ifits_init`:
  Initialise the interface for file reading.
Return:
  Address of the interface.
******************************************************************************/
IFFILE *ifits_init(void);

/******************************************************************************
Function `ifits_destroy`:
  Deconstruct the interface for file reading.
Arguments:
  * `ifile`:    structure for file reading.
******************************************************************************/
void ifits_destroy(IFFILE *ifile);

/******************************************************************************
Function `ifits_newfiles`:
  Open a set of new files for reading.
Arguments:
  * `ifile`:    interface for file reading;
  * `fnames`:   names of the file to be read from;
  * `ninput`:   number of input files to be read from.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int ifits_newfiles(IFFILE *ifile, const char **fnames, const int ninput);

/******************************************************************************
Function `ifits_readlines`:
  Read multiple records (lines) from the input file.
Arguments:
  * `ifile`:    interface for file reading;
  * `nline`:    number of lines to be read at once.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int ifits_readlines(IFFILE *ifile, const size_t nline);

#endif          /* WITH_CFITSIO */

#endif
