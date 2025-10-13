/**
 *
 * @file test_ddot.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Functions to test the ddot variants on lapack format.
 *
 * @version 0.2.0
 * @author Mathieu Faverge
 * @date 2021-09-23
 *
 */
#include "algonum_int.h"
#include <stdio.h>
#include <stdlib.h>


/**
 * Ajoute une ligne de résultats (valeurs numériques) au fichier, séparées par des espaces.
 * 
 * @param filename Nom du fichier où écrire
 * @param values Tableau de double contenant les mesures
 * @param N Nombre de valeurs
 * 
 * @return 0 si succès, -1 en cas d’erreur
 */
int append_result_line(const char *filename, const double *values, size_t N) {
    FILE *file = fopen(filename, "a");  // Ouverture en mode ajout
    if (!file) {
        perror("Erreur lors de l'ouverture du fichier");
        return -1;
    }

    for (size_t i = 0; i < N; ++i) {
        fprintf(file, "%.6f", values[i]);  // Affiche avec 6 décimales
        if (i < N - 1) {
            fprintf(file, " ");  // Séparateur : espace
        }
    }
    fprintf(file, "\n");
    fclose(file);
    return 0;
}

double flops_ddot(int N){
  return flops_dgemm( 1, 1, N );
}

int
testone_ddot( ddot_fct_t ddot, int N, int check )
{
    int     rc = 0;
    double *A, *B, *C;
    int     lda, ldb, ldc;
    double  alpha = 1.0;
    double  beta = 0.0;
    int     seedA = random();
    int     seedB = random();
    int     seedC = random();
    const char *filename = "resultats.txt";
    double mesures[] = {0,0};
    perf_t  start, stop;

    double gflops;
    double flops = flops_ddot( N );
    lda = 1;
    ldb = 1;
    ldc = 1;

    /* Allocate A, B, and C */
    A = malloc( lda * N * sizeof(double) );
    B = malloc( ldb * N * sizeof(double) );
    C = malloc( ldc * 1 * sizeof(double) );

    /* Fill the matrices with random values */
    CORE_dplrnt( 0, 1, N, A, lda, 1, 0, 0, seedA );
    CORE_dplrnt( 0, 1, N, B, ldb, N, 0, 0, seedB );
    CORE_dplrnt( 0, 1, 1, C, ldc, 1, 0, 0, seedC );

    /* Calculate the product */
    perf( &start );
    *C = ddot( N, A, lda, B, ldb );
    perf( &stop );

    perf_diff( &start, &stop );
    if ( flops > 0. ) {
        gflops = perf_gflops( &stop, flops );
    }
    else {
        gflops = 0.;
    }

    /* Check the solution */
    if ( check ) {
        double *Cinit = malloc( ldc * 1  * sizeof(double) );
        CORE_dplrnt( 0, 1, 1, Cinit, ldc, 1, 0, 0, seedC );

        rc = check_ddot( N, A, lda, B, ldb, Cinit, C );

        if ( rc == ALGONUM_SUCCESS) {
            printf( "N= %4d: %le GFlop/s: " ALGONUM_COLOR_GREEN
                    "SUCCESS\n" ALGONUM_COLOR_RESET, N, gflops);
            mesures[0] = N;
    	    mesures[1] = gflops;
            append_result_line(filename, mesures, 2);
        }	
        else {
            printf( "N= %4d: %le GFlop/s: " ALGONUM_COLOR_RED
                    "FAIL\n" ALGONUM_COLOR_RESET, N, gflops);
        }
        free( Cinit );
    }
    else {
	mesures[0] = N;
    	mesures[1] = gflops;
    	append_result_line(filename, mesures, 2);
        printf( "N= %4d: %le GFlop/s\n", N, gflops );
    }
  
    free( A );
    free( B );
    free( C );

    return rc;
}

/**
 * @brief Function to test a series of tests on a ddot function.
 *
 * @param[in] ddot
 *          The function pointer to a ddot operation on lapack format to test.
 *
 * @retval 0, on success
 * @retval The number of failures otherwise.
 *
 */
int
testall_ddot( ddot_fct_t ddot )
{
    int all_N[] = { 0, 3, 5, 17, 64 };

    int nb_N = sizeof( all_N ) / sizeof( int );

    int im, in, ik, m, n, k;
    int nbfailed = 0;
    int nbpassed = 0;
    int nbtests = nb_N;
    CBLAS_TRANSPOSE tA, tB;

    for( in = 0; in < nb_N; in ++ ) {
        n = all_N[in];
        if ( testone_ddot(ddot, n, 1) != ALGONUM_SUCCESS ) {
            nbfailed += 1;
        }
        nbpassed++;
        fprintf( stdout, "Test %4d / %4d\n", nbpassed, nbtests );
    }

    if ( nbfailed > 0 ) {
        printf( ALGONUM_COLOR_RED "\n %4d tests failed out of %d\n" ALGONUM_COLOR_RESET,
                nbfailed, nbtests );
    }
    else {
        printf( ALGONUM_COLOR_GREEN "\n Congratulations all %4d tests succeeded\n" ALGONUM_COLOR_RESET,
                nbtests );
    }
    return nbfailed;
}

int
testwarm_ddot( ddot_fct_t ddot, int N )
{
    int     rc = 0;
    double *A, *B, C;
    int     lda, ldb, ldc, i;
    double  alpha = 1.0;
    double  beta = 0.0;
    int     seedA = random();
    int     seedB = random();
    int64_t nbiter;
    perf_t  start, stop;

    double gflops;
    double flops = flops_ddot( N );
    lda = 1;
    ldb = 1;
    ldc = 1;

    /* Allocate A, B, and C */
    A = malloc( lda * N * sizeof(double) );
    B = malloc( ldb * N * sizeof(double) );

    /* Fill the matrices with random values */
    CORE_dplrnt( 0, 1, N, A, lda, 1, 0, 0, seedA );
    CORE_dplrnt( 0, 1, N, B, ldb, N, 0, 0, seedB );

    /* Estimate the number of iteration to get a 2 seconds run @ 10GFLop/s */
    nbiter = ( 2. / flops ) * 1e10;
    nbiter = my_imin( 1000000, nbiter );
    nbiter = my_imax( 3, nbiter );

    /* Calculate the product */
    perf( &start );
    for( i=0; i<nbiter; i++ ) {
        C += ddot( N, A, lda, B, ldb );
    }
    perf( &stop );

    perf_diff( &start, &stop );
    if ( flops > 0. ) {
        gflops = perf_gflops( &stop, flops * nbiter );
    }
    else {
        gflops = 0.;
    }

    printf( "N= %4d : %le GFlop/s\n (%7ld iterations)", N, gflops, nbiter );

    free( A );
    free( B );

    return rc;
}
