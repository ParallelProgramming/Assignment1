/*
 ============================================================================
 Name        : MaxPrimeGap.cpp
 Author      : lovejoy (Tianran Wang)
 Version     : 1.8
 Description : Finds the largest gap between consecutive primes using
               parallel approach
 ============================================================================
 */
#include <iostream>
#include <mpi.h>
#include <gmp.h>
#include <stdlib.h>

using namespace std;

// check if the result is correct or not -> https://en.wikipedia.org/wiki/Prime_gap

long long mpz_get_ll(mpz_t t);
void mpz_set_ll(mpz_t t, long long l);

int main(int args, char **argv)
{
    int pro_size; // size of processes
    int pro_rank; // rank of current process
    double start_time, end_time;
    MPI_Status status;
    long long upper_limit;

    if (args > 1)
    {
        upper_limit = atoll(argv[1]);
    } else
    {
        upper_limit = 1000000000;
    }

    if (MPI_Init(&args, &argv) != MPI_SUCCESS)
    {
        cout << "mpi init error" << endl;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &pro_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &pro_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    long long quotient = upper_limit / pro_size; // e.g. q=3 u=20 p=6
    int remainder = pro_size + upper_limit % pro_size; // e.g. r=6+(-4)
    int shift_sum = pro_rank < remainder ? pro_rank : remainder;
    int shift = pro_rank < remainder ? 1 : 0;
    long long start = pro_rank * quotient + shift_sum + 1; // first value of current process
    long long next_start = start + quotient + shift; // first value of next process if it exists

    long long curr, next_prime, pre_prime, left_prime, right_prime, max_gap, gap;
    mpz_t gmp_curr, gmp_next_prime;
    mpz_init(gmp_curr);
    mpz_init(gmp_next_prime);
    for (curr = start, pre_prime = 0; curr < next_start && curr != pre_prime; curr = next_prime)
    {
        if (curr !=start) // each process will calc the gap of the end,
        {
            pre_prime = next_prime; // ignore the gap at the beginning
        }

        mpz_set_ll(gmp_curr, curr);
        mpz_nextprime(gmp_next_prime, gmp_curr); // gmp build-in method
        next_prime = mpz_get_ll(gmp_next_prime);

        if (next_prime > upper_limit) // if it exceed the upper limit
        {
            next_prime = pre_prime;
        }

        if (curr !=start) // calc gap for each 2 prime numbers
        {
            gap = next_prime - pre_prime;
        }

        if (gap > max_gap) // update max gap of each process
        {
            max_gap = gap;
            left_prime = pre_prime;
            right_prime = next_prime;
            gap = 0;
        }
    }

    if (pro_rank == 0)
    {
        long long global_max_gap = max_gap; // get results from itself
        long long global_left_prime = left_prime;
        long long global_right_prime = right_prime;

        for (int i = 1; i < pro_size; i++) // update them with data from other processes
        {
            MPI_Recv(&max_gap, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&left_prime, 1, MPI_LONG_LONG, i, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&right_prime, 1, MPI_LONG_LONG, i, 2, MPI_COMM_WORLD, &status);

            if (max_gap > global_max_gap)
            {
                global_max_gap = max_gap;
                global_left_prime = left_prime;
                global_right_prime = right_prime;
            }
        }

        std::cout << "upper limit is " << upper_limit
                  << "\nmax prime gap is " << global_max_gap
                  << "\nleft prime is " << global_left_prime
                  << "\nright prime is " << global_right_prime << std::endl;
    } else
    {
        MPI_Send(&max_gap, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&left_prime, 1, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
        MPI_Send(&right_prime, 1, MPI_UNSIGNED_LONG, 0, 2, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    if (pro_rank == 0)
    {
        std::cout << "runtime is " << end_time - start_time << "s" << std::endl;
    }

    MPI_Finalize();
    return 0;
}

long long mpz_get_ll(mpz_t t)
{
    long long val = 0;
    mpz_export(&val, 0, -1, sizeof val, 0, 0, t);
    return val;
}

void mpz_set_ll(mpz_t t, long long l)
{
    mpz_import(t, 1, -1, sizeof l, 0, 0, &l);
}