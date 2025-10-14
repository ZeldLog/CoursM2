/**
 *
 * @file check_ddot.c
 *
 * @copyright 2019-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Binary checking the numerical validity of the ddot function
 *
 * @version 0.2.0
 * @author Mathieu Faverge
 * @date 2021-09-30
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include "algonum.h"

int init_result_file(const char *filename, const char **var_names, size_t N) {
    FILE *file = fopen(filename, "w");  //  ^icriture (efface le contenu pr  c  dent)
    if (!file) {
        perror("Erreur lors de l'ouverture du fichier");
        return -1;
    }

    for (size_t i = 0; i < N; ++i) {
        fprintf(file, "%s", var_names[i]);
        if (i < N - 1) {
            fprintf(file, " ");  // S  parateur : espace
        }
    }
    fprintf(file, "\n");
    fclose(file);
    return 0;
}

int main( int argc, char **argv )
{
    option_t *options = &global_options;
	
    const char *var_names[] = {"taille", "GFlop/s"};
    init_result_file("resultats.txt", var_names, 2);	

    algonum_init( argc, argv, options, ALGO_DDOT );

    CBLAS_TRANSPOSE tA, tB;
    
    int i = 4;
    while(i < 5) {
        if ( testone_ddot(options->fct->fctptr, i, 1) != ALGONUM_SUCCESS ) {
            printf("fail\n");
        }
	i++;
        
    }

    algonum_exit( options, ALGO_DDOT );

    return EXIT_SUCCESS;
}
