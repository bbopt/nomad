// Blackbox evaluation for a bbob-constrained pb

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <time.h>
#include <assert.h>


#include "coco.h"


#define max(a,b) ((a) > (b) ? (a) : (b))

const size_t DIMENSION = 10;
const size_t INSTANCE = 1;
const size_t FUNCTION_INDEX = 1;
const char COCO_SUITE_NAME[] = "bbob-constrained";

const size_t NB_CONS_FOR_PARAMS = 1;

/**
 * The main method calls the example experiment of the bbob constrained suite.
 */
int main(int argc, char **argv)
{
    
    
    double x[DIMENSION];
    if (argc >= 2)
    {
        FILE * fp;
        char * line;
        size_t len = 0;
        fp = fopen(argv[1], "r");
        if (fp == NULL)
        {
            printf("Input file cannot be opened");
            return -1;
        }
        for (size_t i = 0; i < DIMENSION; i++ )
        {
            fscanf(fp, "%lf", &x[i]);
            if (feof(fp))
            {
                printf("Error reading double");
                return -1;
            }
            // printf("%lf ",x[i]);
        }
        // printf("\n");
        fclose(fp);
    }
    else
    {
        printf("Input file required");
        return -1;
    }
    
    coco_suite_t *suite = coco_suite(COCO_SUITE_NAME, "", "");
    
    coco_problem_t *PROBLEM = coco_suite_get_problem_by_function_dimension_instance(suite,FUNCTION_INDEX,DIMENSION,INSTANCE);

    
    const size_t nb_obj = coco_problem_get_number_of_objectives(PROBLEM);
    if (nb_obj != 1)
    {
        printf("Single objective only");
        return -1;
    }
    const size_t nb_cons = coco_problem_get_number_of_constraints(PROBLEM);
    
    if (nb_cons != NB_CONS_FOR_PARAMS)
    {
        printf("Inconsistent number of constraints for this problem.");
        return -1;
    }
    
    
    // Ready for evaluation
    double f;
    double g[nb_cons];
    
    coco_evaluate_function(PROBLEM, x, &f);
    printf("%f ",f);
    if ( nb_cons >=1 )
    {
        coco_evaluate_constraint(PROBLEM, x, g);
        for (size_t i=0 ; i < nb_cons ; i++)
        {
            printf("%f ",g[i]);
        }
        printf("\n");
    }
    
    coco_suite_free(suite);
    
    return 0;
}

