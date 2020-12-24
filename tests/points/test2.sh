#!/bin/sh

./h2p-ricpad \
    fixed \
    --Dmin 5 \
    --Dmax 40 \
    --d 1 \
    --U0 -0.04 --R0 10  --A0 -41 \
    --m 3 --s 1 \
    --ndigits 600 \
    --tol "1E-90" \
    --h "1E-200" \
    --use-E 

