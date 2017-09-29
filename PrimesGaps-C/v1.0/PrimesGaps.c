/*
 ============================================================================
 Name        : PrimesGaps.c
 Author      :
 Version     : 1.0
 Description : Finds the largest gap between consecutive primes using
               parallel approach
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <mpi.h>
#include <math.h>
#include <gmp.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

const unsigned long long RANGE = 1000000000;
const int MASTER = 0;

int max_gap(int rank, int p_num);
void reduce_gaps();
void p_printf(const char* msg, ...);
unsigned long long mpz_get_ull(mpz_t t);
void mpz_set_ull(mpz_t t, unsigned long long l);

int p_num, rank;
double start_time, run_time;
time_t mytime;

int gap;
unsigned long long p_range, low_prime;

int main(int argc,char**argv) {
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	p_range = floor(RANGE/p_num);

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == MASTER){
		start_time = MPI_Wtime();
		p_printf("Starting");
	}

	gap = max_gap(rank, p_num);

	reduce_gaps();

	if(rank == MASTER){
		p_printf("Largest gap found: %d, between %lli and %lli.", gap, low_prime, low_prime + gap);

		 run_time = MPI_Wtime();
		 p_printf("Run time was %f seconds",run_time-start_time);
		 p_printf("End");
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

int max_gap(int rank, int p_num)
{
	// Boundaries
	mpz_t start;
	mpz_t end;
	mpz_init(start);
	mpz_set_ull(start, rank*p_range + MIN(RANGE % p_num, rank));
	mpz_init(end);

	// adjust self range
	if (RANGE % p_num > rank){
		p_range += 1;
	}
	mpz_add_ui (end, start, p_range);

	p_printf("Working on range: %lli to %lli", mpz_get_ull(start), mpz_get_ull(end));

	mpz_t curr_prime, next_prime, gap, max_gap;
	mpz_init(curr_prime);
	mpz_init(next_prime);
	mpz_init(gap);
	mpz_init(max_gap);

	mpz_nextprime(curr_prime, start);

	while(mpz_cmp(end, curr_prime) > 0){
		mpz_nextprime(next_prime, curr_prime);
		mpz_sub(gap, next_prime, curr_prime);

		if (mpz_cmp(gap, max_gap) > 0 && mpz_cmp(next_prime, end) <= 0){
			mpz_set(max_gap, gap);
			low_prime = mpz_get_ull(curr_prime);
		}

		mpz_set(curr_prime, next_prime);
	}

	return mpz_get_ui(max_gap);
}

void reduce_gaps(){
	unsigned long long data[2];

	if (rank != MASTER)
	{
		data[0] = low_prime;
		data[1] = (unsigned long long)gap;
		p_printf("Sending: max gap - %d, low prime - %lli.",  gap, low_prime);
		MPI_Send(&data, 2, MPI_LONG_LONG, MASTER, 0, MPI_COMM_WORLD);
	}
	else{
		int i;
		for (i = 1; i < p_num; i++){
			MPI_Recv(&data, 2, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			p_printf("Received from process %d: max gap - %d, low prime - %lli.", i,  data[1], data[0]);

			if( data[1] > gap){
				gap =  data[1];
				low_prime =  data[0];
			}
		}
	}
}

void p_printf(const char* format, ...){
	mytime = time(NULL);
	printf("%s Process %d / %d: ", strtok(ctime(&mytime), "\n"), rank, p_num);

	va_list arg;
	va_start(arg, format);
	vfprintf (stdout, format, arg);
	va_end(arg);
	printf("\n");
}

unsigned long long mpz_get_ull(mpz_t t)
{
    unsigned long long val = 0;
    mpz_export(&val, 0, -1, sizeof val, 0, 0, t);
    return val;
}

void mpz_set_ull(mpz_t t, unsigned long long l)
{
    mpz_import(t, 1, -1, sizeof l, 0, 0, &l);
}

