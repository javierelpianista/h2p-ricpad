#!/bin/sh

# This script allows to compute the equilibrium internuclear distance, 
# separation constant and electronic + nuclear energy with D = 2..15.
#
# No output to any file is produced.
#
# h2p-ricpad should be linked or copied over to this folder

./h2p-ricpad \
    minimum \
    --Dmin 2 \
    --Dmax 15 \
    --d 0 \
    --U0 -0.6 \
    --R0 2.0  \
    --A0 0.81 \
    --m 0 \
    --s 0 \
    --ndigits 1200 \
    --tol "1E-90" \
    --h "1E-200" \
    --hd "1E-400" \
    --no-log
