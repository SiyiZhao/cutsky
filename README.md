# cutsky

![GitHub](https://img.shields.io/github/license/cheng-zhao/cutsky.svg)
![Codacy grade](https://img.shields.io/codacy/grade/b8e918a70b3b486e982cf55fb6104b88.svg)

## Table of Contents

- [Introduction](#introduction)
- [Compilation](#compilation)
- [Running](#running)
- [References](#references)

## Introduction

This toolkit constructs survey-like cut-sky galaxy catalogues from simulations in periodic cubic boxes. The workflow involves:
-   Tiling the simulation volume by duplicating the box with periodic boundary conditions, to enclose the full survey volume;
-   Converting comoving coordinates to observational coordinates (right ascension, declination, redshift) assuming the observer is located at the origin (0, 0, 0);
-   Applying redshift-space distortions using peculiar velocities;
-   Trimming to the specified survey footprint;
-   Applying a radial selection to match a target number density as a function of redshift.

The input catalogue may be either an ASCII file with the first 6 columns representing comoving coordinates and peculiar velocities of tracers (`x`, `y`, `z`, `vx`, `vy`, `vz`), or a FITS table with equivalent columns. The columns of the output cut-sky catalogues are listed as follows:

| Column        | Description                                                                                                |
|---------------|------------------------------------------------------------------------------------------------------------|
| `RA`          | Right Ascension                                                                                            |
| `DEC`         | Declination                                                                                                |
| `Z`           | Observed redshift                                                                                          |
| `Z_COSMO`     | &ldquo;Real-space&rdquo; redshift (without peculiar velocity effects)                                      |
| `STATUS`      | *(Optional)* Bitmask for extra selections: `1` = within an extra footprint; `2` = kept in radial selection |
| `NZ`          | *(Optional)* Comoving number density                                                                       |
| `RAN_NUM_0_1` | *(Optional)* Random number in [0, 1) for radial selection                                                  |

The main difference between this code and [`make_survey`](https://github.com/mockFactory/make_survey) is that it directly duplicates the cubic simulation box, rather than remapping it to a cuboid using [`BoxRemap`](http://mwhite.berkeley.edu/BoxRemap/). This approach typically requires smaller box sizes to fill a survey volume and readily enables the construction of cut-sky catalogues from small simulations by reusing the same box. In addition, `cutsky` is highly optimised and parallelised, making it suitable for the production of a large number of mock catalogues. It has been used for generating the DESI DR1 EZmocks<sup>[\[1\]](#ref1)</sup>.

This program is compliant with the ISO C99 and IEEE POSIX.1-2008 standards, and no external library is mandatory. Thus it is compatible with most modern C compilers and operating systems. Optional support for the `FITS` file format is available via the [`cfitsio` library](https://github.com/HEASARC/cfitsio) (see also [Compilation](#compilation)).

This software is written by Cheng Zhao (&#36213;&#25104;), and is distributed under the [MIT license](LICENSE.txt). If you use this program in research work that results in publications, please cite [\[1\]](#ref1).

## Compilation

The build process of `cutsky` is based on the `make` utility. Compilation options can be customised in the [`options.mk`](options.mk) file, including the flag to enable FITS-format support.

If FITS support is enabled, `cutsky` requires the [`cfitsio` library](https://github.com/HEASARC/cfitsio). If it is installed in a non-standard location, please specify the path via the `CFITSIO_DIR` entry in [`options.mk`](options.mk#L10). The compiler will then look for `fitsio.h` in `CFITSIO_DIR/include`, and the library file (`libcfitsio.so`, `libcfitsio.a`, or `libcfitsio.dylib`) in `CFITSIO_DIR/lib`.

Once configured, compile the program with:

```bash
make
```

The following command clear the compiled executable:
```bash
make clean
```

## Running

The default executable `CUTSKY` can be run with:

```bash
./CUTSKY
```

By default, it reads the configuration file `cutsky.conf` in the current working directory. A template configuration file can be generated with:

```bash
./CUTSKY -t > cutsky.conf
```

Alternatively, configuration options may be passed via command-line options, including the location of a custom configuration file. These override any values set in the configuration file. A summary of available command-line options can be displayed via

```bash
./CUTSKY -h
```

A detailed explanation of all configuration parameters is provided in [`cutsky.conf`](cutsky.conf).

The scripts for generating survey footprints in Mangle polygon format from DESI tiles are provided in the [`scripts`](scripts) directory.

## References

<span id="ref1">\[1\]</span> Zhao et al., 2025, Mock catalogues with survey realism for the DESI DR1, *in preparation*.
