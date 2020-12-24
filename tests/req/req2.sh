#!/bin/sh

./h2p-ricpad minimum \
    --Dmin 2 \
    --U0 -0.602 \
    --A0 0.8 \
    --R0 2.0 \
    --d 2 \
    --ndigits 1000 \
    --h "1E-200" \
    --hd "1E-400" \
    --tol "1E-90" \
    --nr-max-iter 10000 \
    --output-file d2.out \
    --log-nr 
