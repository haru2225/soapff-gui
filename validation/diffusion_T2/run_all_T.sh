#!/bin/sh
# 300, 330 and 366 K one after the other (np 3 each: three concurrent jobs are 5x slower per job on this laptop)
cd "$(dirname "$0")"
for T in 300 330 366; do NP=3 ./make_T.sh $T > T$T.log 2>&1; done
