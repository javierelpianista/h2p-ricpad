#!/bin/sh

./h2p-ricpad minimum \
    --Dmin 2 \
    --U0 -0.602 \
    --A0 0.8 \
    --R0 2.0 \
    --d 2 \
    --ndigits 2000 \
    --h "1E-400" \
    --hd "1E-800" \
    --tol "1E-120" \
    --nr-max-iter 10000 \
    --output-file d2.out \
    --log-nr 
