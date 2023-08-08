#include "numerics_disease.h"

FILE* initialize_file(char* filename, void* params){
    double* param = (double*)params; //param contains signature parameter which are included in the filenames making each file identifiable 
    int signature = (int)param[0];
    double p = param[1];
    char fn[50];
    sprintf(fn, filename, signature,p);
    FILE* fp = fopen(fn, "w");
    if(fp == NULL){
        printf("File could not be generated! \n");
        return NULL;
    }
    else{
        return fp;
    }
}

double mean(double* array, int len){
    double sum = 0;
    for(int i = 0; i < len; i++){
        sum += array[i];
    }
    return sum / len;
}

double standard_deviation(double* array, int len){
	double Mean = mean(array, len);
	double residues = 0;
	for(int i = 0; i < len; i++){
		residues += pow((array[i] - Mean), 2);
	}
	double res = (double)residues;
	return sqrt(res)/ (double)(len - 1);
}

//Checks whether a pointer has been successfully initialized. Returns True if correctly init. and False if the return value is NULL
int check_pointer(void* pointer){
	if(pointer == NULL){ return 0;}
	else{ return 1;}
}

double ran_num(gsl_rng* generator){
	return gsl_rng_uniform(generator);
}
//generates a random float between min and max
int ran_int(gsl_rng* generator, int min, int max){
	return gsl_rng_uniform_int(generator, max - min +1) +min;
}


// An array of randomly generated numbers can be avaluated using the crosscorrelation. The goodness of randomness can then be derived from the results. The more arbitrary the autocorrelation (abscence of a pattern) the more random the generated numbers
void auto_correlation(double* signal, int length){
	FILE* fp = fopen("data/autocorr.csv", "w");
	double* autocorr = (double*) malloc(sizeof(double)*length);
	if(!check_pointer(fp)){printf("File could not be created! \n"); return;}
	if(!check_pointer(autocorr)){ printf("System ran out of memory! \n"); return;}
	
	double signal_mean = mean(signal, length); // Calculating the mean of the signal using the mean method from numerics.c
	for(int k = 0; k < length; k++){
		autocorr[k] = 0.0;
		for(int i = 0; i < length-k; i++){
			autocorr[k] += (signal[i]-signal_mean) * (signal[i+k] - signal_mean);
		}
		autocorr[k] /= length;
		fprintf(fp,"%g, %g\n", signal[k], autocorr[k] );
	}
	free(autocorr);
}

//3D-Plotting of randomly generated triplets. If no hyperplanes can be identified the randomness of generation is sufficiently good
void goodness_of_random(double* signal, int length){
    FILE* fp = fopen("data/goodness_random_triple.csv", "w");
    if(!check_pointer(fp)){printf("File could not be created! \n"); return;} 
    
    int check = 0;
    int len = length;
    while(!(check = (len % 3 == 0))){
        len--;
    }
    for(int i = 0; i < (len /3); i+=3){
        fprintf(fp, "%g, %g, %g\n", signal[i],signal[i+1], signal[i+2]);
    }
    fclose(fp);
}


