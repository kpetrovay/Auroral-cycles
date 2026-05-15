#!/bin/sh
rm -f mcsim_results.dat
for F in 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0
#for F in 3.0 3.5 
do
for M in 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 
#for M in 4.5 5.0 
do
echo "Now running aursim.py with args" `expr $F` `expr $M` 
./aursim.py $F $M 
done
done
