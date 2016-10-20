#!/bin/sh
cd lp_solve_5.5/lpsolve55
sh ccc
cd ../..
cp lp_solve_5.5/lpsolve55/bin/*/* .

cd samtools-0.1.19/
make
cd ../

make
