include options.mk

LIBS = -lm

# Settings for CFITSIO
ifeq ($(strip $(WITH_FITS)), T)
  CFLAGS += -DWITH_CFITSIO
  LIBS += -lcfitsio
  ifneq ($(strip $(CFITSIO_DIR)),)
    LIBS += -L$(strip $(CFITSIO_DIR))/lib
    INCL += -I$(strip $(CFITSIO_DIR))/include
  endif
endif

# Settings for OpenMP
ifeq ($(strip $(USE_OMP)), T)
  LIBS += -DOMP -fopenmp
endif

INCL += -Isrc -Iio -Ilib -Imath -Iprand/src/header
SRCS = $(wildcard src/*.c io/*.c lib/*.c math/*.c prand/src/*.c)
EXEC = CUTSKY

all:
	$(CC) $(CFLAGS) -o $(EXEC) $(SRCS) $(LIBS) $(INCL)

clean:
	rm $(EXEC)
