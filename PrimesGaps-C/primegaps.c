/*
 ============================================================================
 Name        : PrimeGaps.c
 Author      :
 Version     : 1.1
 Description : Finds the largest gap between consecutive primes using
               parallel approach
 ============================================================================
 */

#include "primegaps.h"

int main(int argc,char**argv) {
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if( argc == 2 ) {
		mpz_set_ull(range, atoll(argv[1]));
	}
	else {
		 mpz_set_ull(range, DEFAULT_RANGE); /* default is 10^9 */
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == MASTER_RANK){
		start_time = MPI_Wtime();
		p_printf("Starting");
	}

	setup();  /* calculate the work-range */

	max_gap();  /* find max gap */

	reduce_gaps();  /* send all gaps to master and find the largest */

	if(rank == MASTER_RANK){
		p_printf("Largest gap found: %d, between %lli and %lli.", p_max_gap,
				p_prime, p_prime + p_max_gap);

		 end_time = MPI_Wtime();
		 p_printf("Run time was %f seconds", end_time-start_time);
		 p_printf("End");
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

/*
 * Function: setup
 * ----------------------------
 *   Calculate the range (start and end values) for the current process to work on
 *   Using: n(p) 	 = floor(n/P) ( +1 if p < mod(n,P)) = p_range
 *   		start(p) = p*floor(n/P) + min(p, mod(n,P)) 	= p_start
 *   		where p is rank
 */
void setup(){
	mpz_t rem;

	mpz_init(p_start);
	mpz_init(p_end);
	mpz_init(p_range);
	mpz_init(rem);

	mpz_fdiv_q_ui(p_range, range, p_num); 	/* p_range = floor(n/P) */
	mpz_fdiv_r_ui(rem, range, p_num);     	/* rem = mod(n,P) */

	mpz_mul_ui(p_start, p_range, rank); 	/* p_start = p*floor(n/P) */

	if(mpz_cmp_ui(rem, rank) > 0){ 			/* if mod(n,P) > p */
		mpz_add_ui(p_start, p_start, rank); /* p_start += p */
		mpz_add_ui(p_range, p_range, 1); 	/* increase p_range by 1 */
	}
	else{
		mpz_add(p_start, p_start, rem);  	/* p_start += mod(n,P) */
	}

	mpz_add(p_end, p_start, p_range); 		/* end(p) = start(p) + n(p) */

	p_printf("Working on range: %lli to %lli", mpz_get_ull(p_start), mpz_get_ull(p_end));
}

/*
 * Function: max_gap
 * ----------------------------
 *   Finds and sets p_max_gap and p_prime for the process work-range
 */
void max_gap()
{
	mpz_t curr_prime, next_prime, gap, max_gap;
	mpz_init(curr_prime);
	mpz_init(next_prime);
	mpz_init(gap);
	mpz_init(max_gap);

	/* start from the first prime in the process work-range */
	mpz_nextprime(curr_prime, p_start);

	while(mpz_cmp(p_end, curr_prime) > 0){
		mpz_nextprime(next_prime, curr_prime);
		mpz_sub(gap, next_prime, curr_prime);

		/* if the gap is the largest so far and it's not outside our entire range */
		if (mpz_cmp(gap, max_gap) > 0 && mpz_cmp(next_prime, range) <= 0){
			mpz_set(max_gap, gap); /* update the max */
			p_prime = mpz_get_ull(curr_prime);
		}

		mpz_set(curr_prime, next_prime); /* move on to the next prime */
	}

	p_max_gap = mpz_get_ull(max_gap); /* set the max gap found */
}

/*
 * Function: reduce_gaps
 * ----------------------------
 *   Sends all max gaps found by all non-MASTER processes to MASTER.
 *   Stores the largest gap and it's first occurrence in the MASTER process.
 */
void reduce_gaps(){
	unsigned long long data[2];

	if (rank != MASTER_RANK)
	{
		data[0] = p_max_gap;
		data[1] = p_prime;
		MPI_Send(&data, 2, MPI_LONG_LONG, MASTER_RANK, 0, MPI_COMM_WORLD);
	}
	else{
		p_printf("Found gap - %d, low prime - %lli.", p_max_gap, p_prime);

		MPI_Status status;
		int i;
		for (i = 0; i < p_num - 1; i++){
			MPI_Recv(&data, 2, MPI_LONG_LONG, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			p_printf("Received from process %d: max gap - %d, low prime - %lli.",
					status.MPI_SOURCE, data[0], data[1]);
			if( data[0] > p_max_gap){
				p_max_gap =  data[0]; /* update max gap */
				p_prime =  data[1];
			}
		}
	}
}

/*
 * Function: p_printf
 * ----------------------------
 *   Similar to printf, adds the process identifier and a timestamp.
 */
void p_printf(const char* format, ...){
	t = time(NULL);
	printf("%s Process %d / %d: ", strtok(ctime(&t), "\n"), rank, p_num);

	va_list arg;
	va_start(arg, format);
	vfprintf (stdout, format, arg);
	va_end(arg);
	printf("\n");
}

/*
 * Function: mpz_get_ull
 * ----------------------------
 *   Converts the value of a mpz_t variable to unsigned long long
 *   t: the mpz_t variable
 *
 *   returns: the value of t as unsigned long long
 */
unsigned long long mpz_get_ull(mpz_t t)
{
    unsigned long long val = 0;
    mpz_export(&val, 0, -1, sizeof val, 0, 0, t);
    return val;
}

/*
 * Function: mpz_set_ull
 * ----------------------------
 *   Sets the value of a mpz_t variable from unsigned long long variable
 *   t: the mpz_t variable to set
 *   l: the value to set in t
 */
void mpz_set_ull(mpz_t t, unsigned long long l)
{
    mpz_import(t, 1, -1, sizeof l, 0, 0, &l);
}
