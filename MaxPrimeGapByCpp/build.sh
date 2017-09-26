mpic++ ./MaxPrimeGap.cpp -o MaxPrimeGap -lgmp
mpirun -n 8 -host orc-dev3 ./MaxPrimeGap 1000000000
