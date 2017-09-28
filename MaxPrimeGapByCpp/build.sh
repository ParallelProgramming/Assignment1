# to build
mpic++ ./MaxPrimeGap.cpp -o MaxPrimeGap -lgmp

# to run on local environment
mpirun -np 8 ./MaxPrimeGap 1000000000

# to run on orca cluster
sqsub -q mpi -r 3d -o ofile_9_8.log  -n 8 ./MaxPrimeGap 1000000000
