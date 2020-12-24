#!/bin/sh

./h2p-ricpad \
    fixed \
    --Dmin 2 \
    --Dmax 80 \
    --d 0 \
    --U0 -0.6 \
    --R0 2.0  \
    --A0 0.81 \
    --m 0 \
    --s 0 \
    --ndigits 1200 \
    --tol "1E-90" \
    --h "1E-200" \
    --hd "1E-400"
