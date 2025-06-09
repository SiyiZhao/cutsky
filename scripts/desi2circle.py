#!/usr/bin/env python3
#
# desi2circle.py: this file is part of the cutsky program.
#
# cutsky: cutsky catalogue generator.
#
# Github repository:
#       https://github.com/cheng-zhao/cutsky
#
# Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

import numpy as np
from astropy.table import Table
from desimodel.focalplane import get_tile_radius_deg    # 1.628032452048558287316382120479829609394
import sys

samples = ['Y5', 'Y1_dark', 'Y1_bright']
ftiles = ['/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-main.ecsv',
       '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/tiles-DARK.fits',
       '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/tiles-BRIGHT.fits']

if len(sys.argv) != 2:
  print(f'Usage: {sys.argv[0]} SAMPLE')
  print('  Possible SAMPLE options:', ' '.join(samples))
  exit(1)

try:
  idx = samples.index(sys.argv[1])
except:
  print(f'Unknown SAMPLE: {sys.argv[1]}')
  print('  Possible SAMPLE options:', ' '.join(samples))
  exit(1)

tiles = Table.read(ftiles[idx])
if idx == 0: tiles = tiles[(tiles['PROGRAM'] != 'BACKUP') & tiles['IN_DESI']]

tiles['RADIUS'] = np.full_like(tiles, get_tile_radius_deg())

tiles.meta['comments'] = ['circle 0 1', 'unit d']

ofile = f'{sys.argv[1]}_circle.dat'
formats = {'RA':'%.10g', 'DEC':'%.10g', 'RADIUS':'%.40g'}

tiles['RA','DEC','RADIUS'].write(ofile, overwrite=True,
    format='ascii.no_header', formats=formats)

print(f'Circle-format mask saved to {ofile}')
