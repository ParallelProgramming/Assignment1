mpic++ ./MaxPrimeGap.cpp -o MaxPrimeGap -lgmp

sqsub -q mpi -r 3d -o ofile_9_8.log  -n 8 ./MaxPrimeGap 1000000000
