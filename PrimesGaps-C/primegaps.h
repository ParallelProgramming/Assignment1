/*
 ============================================================================
 Name        : PrimeGaps.c
 Author      :
 Version     : 1.1
 Description : Finds the largest gap between consecutive primes using
               parallel approach
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <gmp.h>

void setup();
void max_gap();
void reduce_gaps();
void p_printf(const char* msg, ...);
unsigned long long mpz_get_ull(mpz_t t);
void mpz_set_ull(mpz_t t, unsigned long long l);

const int MASTER_RANK 	= 0;
const int DEFAULT_RANGE = 1000000000;

int 				p_num, rank, p_max_gap;
unsigned long long 	p_prime;
mpz_t 				range, p_range, p_start, p_end;
double 				start_time, end_time;
time_t 				t;
