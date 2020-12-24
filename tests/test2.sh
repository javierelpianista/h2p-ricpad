#!/bin/sh

./h2p-ricpad \
    fixed \
    --Dmin 3 \
    --Dmax 80 \
    --d 1 \
    --U0 -0.05389 \
    --R0 10  \
    --A0 -29.37 \
    --m 4 \
    --s 1 \
    --ndigits 600 \
    --tol "1E-90" \
    --h "1E-200" \
    --use-E 

