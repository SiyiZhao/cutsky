#!/bin/bash
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

if [ $# != 2 ]; then
  echo "Usage: $0 INPUT RES"
  exit
fi

IFILE=$1
RES=$2          # 7 for DESI

# Mangle package available at https://space.mit.edu/~molly/mangle/download/
MANGLE_BIN=$HOME/mangle2.2/bin

$MANGLE_BIN/pixelize -vn -p+40 -Ps0,${RES} -ic1d -opd "$IFILE" "${IFILE}.pix${RES}.ply" || exit

$MANGLE_BIN/snap -vn -p+40 -ipd -opd "${IFILE}.pix${RES}.ply" "${IFILE}.snap${RES}.ply" || exit

$MANGLE_BIN/balkanize -p+40 -Bl -ipd -opd "${IFILE}.snap${RES}.ply" "${IFILE}.balk${RES}.ply" || exit

$MANGLE_BIN/unify -p+40 -ipd -opd "${IFILE}.balk${RES}.ply" "${IFILE}_final_res${RES}.ply" || exit

rm "${IFILE}.pix${RES}.ply" "${IFILE}.snap${RES}.ply" "${IFILE}.balk${RES}.ply"
