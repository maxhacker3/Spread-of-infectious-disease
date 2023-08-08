#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

double mean(double*, int);

double standard_deviation(double*, int);

int check_pointer(void*);

double ran_num(gsl_rng*);

int ran_int(gsl_rng*, int, int);

void auto_correlation(double*, int);

void goodness_of_random(double*, int);
    
FILE* initialize_file(char*, void*);
