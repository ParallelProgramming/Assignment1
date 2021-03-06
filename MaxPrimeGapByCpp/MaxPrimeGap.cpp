/*
 ============================================================================
 Name        : MaxPrimeGap.c
 Author      : lovejoy (Tianran Wang)
 Description : To find the largest gap between consecutive primes using
               parallel approach.
 ============================================================================
 */
#include <iostream>
#include <mpi.h>
#include <gmp.h>
#include <stdlib.h>

using namespace std;

// check if the result is correct or not -> https://en.wikipedia.org/wiki/Prime_gap

int main(int args, char **argv)
{
    int pro_size_int; // size of processes
    int pro_rank_int; // rank of current process
    double start_time, end_time;
    MPI_Status status;
    mpz_t upper_limit;

    if (args > 1)
    {
        mpz_init_set_ui(upper_limit, atoll(argv[1]));
    } else
    {
        mpz_init_set_ui(upper_limit, 1000000000);
    }

    if (MPI_Init(&args, &argv) != MPI_SUCCESS)
    {
        cout << "mpi init error" << endl;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &pro_size_int);
    MPI_Comm_rank(MPI_COMM_WORLD, &pro_rank_int);
    start_time = MPI_Wtime();

    mpz_t pro_size; // size of processes
    mpz_init_set_ui(pro_size, pro_size_int);

    mpz_t pro_rank; // rank of current process
    mpz_init_set_ui(pro_rank, pro_rank_int);

    mpz_t quotient, remainder; // q=u/p , r=mod(u,p)-p
    mpz_init(quotient);
    mpz_init(remainder);
    mpz_cdiv_qr(quotient, remainder, upper_limit, pro_size); // e.g. q=3 r=-4 u=20 p=6
    mpz_add(remainder, remainder, pro_size); // e.g. r=-4+6

    mpz_t multiple; // multiple = pro_rank * quotient, for calculating start
    mpz_init(multiple);
    mpz_mul(multiple, pro_rank, quotient);

    mpz_t shift_sum, shift; // shift_sum = min(pro_rank, remainder)
    if (mpz_cmp(pro_rank, remainder) < 0) // pro_rank < remainder
    {
        mpz_init_set(shift_sum, pro_rank);
        mpz_init_set_ui(shift, 1);
    } else // pro_rank >= remainder
    {
        mpz_init_set(shift_sum, remainder);
        mpz_init_set_ui(shift, 0);
    }

    mpz_t start; // first value of current process
    mpz_init(start);
    mpz_add(start, multiple, shift_sum);// start == multiple + shift_sum

    mpz_t next_start; // first value of next process if it exists
    mpz_init(next_start);
    mpz_add(next_start, start, quotient);
    mpz_add(next_start, next_start, shift); // next_start = start + quotient + shift

    mpz_t cur_prime, next_prime, left_prime, right_prime, max_gap, gap;
    mpz_init(cur_prime);
    mpz_init(next_prime);
    mpz_init(left_prime);
    mpz_init(right_prime);
    mpz_init(max_gap);
    mpz_init(gap);

    mpz_nextprime(cur_prime, start);
    while (mpz_cmp(cur_prime, next_start) <= 0)
    {
        mpz_nextprime(next_prime, cur_prime);
        if (mpz_cmp(next_prime, upper_limit) > 0)
        {
            mpz_set_ui(gap, 0);
        } else
        {
            mpz_sub(gap, next_prime, cur_prime);
        }
        if (mpz_cmp(gap, max_gap) > 0) // update max gap of each process
        {
            mpz_set(max_gap, gap);
            mpz_set(left_prime, cur_prime);
            mpz_set(right_prime, next_prime);
            mpz_set_ui(gap, 0);
        }
        mpz_set(cur_prime, next_prime);
    }

    long long int upper_limit_int = mpz_get_ui(upper_limit);
    long long int max_gap_int = mpz_get_ui(max_gap);
    long long int left_prime_int = mpz_get_ui(left_prime);
    long long int right_prime_int = mpz_get_ui(right_prime);

    if (pro_rank_int == 0)
    {
        long long int global_max_gap_int = max_gap_int; // get results from itself
        long long int global_left_prime_int = left_prime_int;
        long long int global_right_prime_int = right_prime_int;

        for (int i = 1; i < pro_size_int; i++) // update them with data from other processes
        {
            MPI_Recv(&max_gap_int, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&left_prime_int, 1, MPI_UNSIGNED_LONG, i, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&right_prime_int, 1, MPI_UNSIGNED_LONG, i, 2, MPI_COMM_WORLD, &status);

            if (max_gap_int > global_max_gap_int)
            {
                global_max_gap_int = max_gap_int;
                global_left_prime_int = left_prime_int;
                global_right_prime_int = right_prime_int;
            }
        }

        std::cout << "upper limit is " << upper_limit_int
                  << "\nmax prime gap is " << global_max_gap_int
                  << "\nleft prime is " << global_left_prime_int
                  << "\nright prime is " << global_right_prime_int << std::endl;
    } else
    {
        MPI_Send(&max_gap_int, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&left_prime_int, 1, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
        MPI_Send(&right_prime_int, 1, MPI_UNSIGNED_LONG, 0, 2, MPI_COMM_WORLD);
    }

    end_time = MPI_Wtime();
    if (pro_rank_int == 0)
    {
        std::cout << "runtime is " << end_time - start_time << "s" << std::endl;
    }

    MPI_Finalize();
    return 0;
}