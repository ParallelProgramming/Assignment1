mpic++ ./MaxPrimeGap.cpp -o MaxPrimeGap -lgmp
mpirun -np 8 ./MaxPrimeGap 1000000000
