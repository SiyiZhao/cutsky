CC = gcc
CFLAGS = -std=c99 -O3 -Wall -flto=auto

USE_OMP = T
WITH_FITS = T  # T for enabling FITS-format outputs

# Directory for CFITSIO (>=4.2) library
# The corresponding header file should be in $(CFITSIO_DIR)/include
# The library file should be in $(CFITSIO_DIR)/lib
CFITSIO_DIR = 
