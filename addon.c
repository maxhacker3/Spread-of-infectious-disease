#include "numerics_disease.h"

typedef enum{
    SUSCEPTIBLE = 1,
    INFECTED = 2,
    RECOVERED = 3,
    VACCINATED = 4
} Status;

typedef struct{
    Status status;
}Grid;

typedef struct{
    double p1;
    double p2;
    double p3;
    double p4;
}Probabilities;

typedef void (*initFunction)(int, Grid**, int, gsl_rng*, void*); //initialization function -> either void initialize or void initialize_single



// In the following certain addons are computed which support the understanding of the results and help evaluating the goodness

// 1. Goodness of random numbers using auto-correlation
// 2. Creating of dataset for Animation -> either conventional initialization of grid (void initialize) or hotspot-initialization (void initialize_hotspots)
// 3. Examination of fluctuation for small grids


void initialize(int L, Grid grid[][L+2], int include_vaccination, gsl_rng* generator, void* params){
    Probabilities* probs = (Probabilities*)params;
    for(int i = 1; i < L+1; i++){
        for(int j = 1; j < L+1; j++){
            double vacc_prob = (include_vaccination) ? ran_num(generator) : 1;
            int random = (vacc_prob <= probs->p4)? 4 : ran_int(generator, 1,3);
            switch(random){
                case 1: grid[i][j].status = SUSCEPTIBLE; break;
                case 2: grid[i][j].status = INFECTED; break;
                case 3: grid[i][j].status = RECOVERED; break;
                case 4: grid[i][j].status = VACCINATED; break;
                default: break;
            }
        }
    }
    
}

//initializing a certain number of hotspots
void initialize_hotspots(int L, Grid grid[][L+2], int include_vaccination, gsl_rng* generator, void* params){
    Probabilities* probs = (Probabilities*)params;
    int number = 30; // number of hotspots
    for(int i = 1; i < L+1; i++){
        for(int j = 1; j < L+1; j++){
            double random_vacc = ran_num(generator);
            grid[i][j].status = (include_vaccination && (random_vacc <= probs->p4))? VACCINATED : SUSCEPTIBLE;
        }
    }
    int num = number;
    for(int n = 0; n < num; n++){
        int random = ran_int(generator,0, pow(L,2)-1);
        //determining the grid correspondent grid indices
        int i = ((int)random /L) + 1;
        int j = (random % L) + 1;
        
        grid[i][j].status = INFECTED;
        
    }
    
}


int count_infected_neighbours(int i1, int j1, int L, Grid grid[][L+2], void* params){
    int i = i1;
    int j = j1;
    int sum = 0;
    sum += (grid[i-1][j].status == INFECTED)? 1: 0;
    sum += (grid[i+1][j].status == INFECTED)? 1: 0;
    sum += (grid[i][j+1].status == INFECTED)? 1: 0;
    sum += (grid[i][j-1].status == INFECTED)? 1: 0;
    
    if(sum != 0){
        return sum;
    }
    else{
        return 0;
    }
}

double calc_probability(double p, int n){
    return 1 - pow(1-p,n); //The probability of infection is calculated according to the binomial distribution. It depends on the number of infected neighbours n
}


void grid_step(double time, int L, Grid grid[][L+2], gsl_rng* generator, void* params){
    Probabilities* probs = (Probabilities*)params;
    double p1 = probs->p1; double p2 = probs->p2; double p3 = probs->p3; double p4 = probs->p4;
    
    // Draws of a random number between 0 and L^2 and converting it to the corresponding grid coordinates.
	int random = ran_int(generator ,0, pow(L,2)-1);
    int i = ((int)random / L) + 1;
    int j = (random % L) + 1;
    
   	double rand_prob = ran_num(generator);// random double between 0 and 1 as probability
    int n; // number of infected neighbours

    switch(grid[i][j].status){ 
    	case 1: //cell susceptible; Number of infected neighbours is determined; probability of infection calculated; cell is updated to infected if applicable
            n = count_infected_neighbours(i,j,L,grid,NULL);
            if(n != 0 && (rand_prob <= calc_probability(p1, n))){
                grid[i][j].status = INFECTED;
            }
            else{;}break;
        case 2://Cell infected; p2 to recover; cell is updated if applicable
            if(rand_prob <= p2){
                grid[i][j].status = RECOVERED;
            }
            else{;}break;
        case 3://cell recovered; p3 to become susceptible; cell is updated if applicable
            if(rand_prob <= p3){
                grid[i][j].status = SUSCEPTIBLE;
            }
            else{;}break;
        case 4://cell is vaccinated and remains vaccinated
            grid[i][j].status = VACCINATED;
        default: break;
    }
}


//the method convieniently prints the dataset for the ANIMATION to a file
void print_to_file(FILE* fp, double time, int L, Grid grid[][L+2]){
    int* flattened = (int*) malloc(sizeof(int) * pow(L,2));
    if(!check_pointer(flattened)) {printf("System ran out of memory! \n"); return;} 
    
    int k = 0;
    for(int i = 1; i < L+1; i++){
        for(int j = 1; j < L+1; j++){
            switch(grid[i][j].status){ //filling the flattened 2d array
                //the printed number equivalents are crucial for correct plotting
        		case 1: flattened[k] = 0;break;
        		case 2: flattened[k] = 2;break;
        		case 3: flattened[k] = 1;break;
        		case 4: flattened[k] = 0;break;
        	}
            k++;
        }
    }
    for(int n = 0; n < pow(L,2)-1; n++){
           fprintf(fp, "%d, ", flattened[n]);
    }
    fprintf(fp, "%d\n", flattened[(int)(pow(L,2)-1)]);
    return;
}


//calculation of the infected cells out of all cells in the grid
double infection_rate(int M, Grid grid[][M+2]){
    int sum = 0;
    for(int i = 1; i < M+1;i++){
        for(int j = 0; j < M+1; j++){
            sum += (grid[i][j].status == INFECTED)? 1 : 0;
        }
    }
    return (double)(sum / pow(M,2));
}

// creating the dataset for the ANIMATION
void animation(int L, double T_max, Grid grid[][L+2], Probabilities probs, int include_vaccination, gsl_rng* generator, initFunction func){
    double param[] = {(double)L, 0}; //signature parameter for the filename
    
    //file management
    FILE* fp_dynamics = initialize_file("data/disease_dynamics.csv", param);
    //next file contains the parameters which are crucial for correct data-processing in Python
    FILE* fp_parameter = fopen("data/animation_parameter.csv", "w");
    if(!check_pointer(fp_dynamics)){printf("File could not be generated!\n"); return;}
    if(!check_pointer(fp_parameter)){printf("File could not be generated!\n"); return;}
    fprintf(fp_parameter, "%f\n%f\n%f\n%f\n", probs.p1, probs.p2, probs.p3, probs.p4);
    fclose(fp_parameter);
    
    
    func(L, (Grid**)grid, include_vaccination, generator, &probs); //initialization
    
    double time = 0;
    double delta_t = 1;
    while(time <= T_max){
        for(int i = 0; i < pow(L,2); i++){
            grid_step(time, L, grid, generator, &probs);
        }
        //print to file
        print_to_file(fp_dynamics, time, L, grid);
        time += delta_t; //incrementing time-step
    }
    fclose(fp_dynamics);
    return;
}

//The function examines the statistical fluctuations for small grids
//It calculates the average infection rate as a function of p1 for N iterations. Eventually, the mean and standard deviation is calculated
void average_noise(int L,int iterations, gsl_rng* generator, void* params){
	Probabilities probs = *(Probabilities*)params;
	    double delta_p = 0.02; // step width of p
	    int N = 1 / 0.02;
	    Grid(*grid)[L+2] = calloc(sizeof(Grid[L+2][L+2]),1);
	    if(!check_pointer(grid)){printf("The grid could not be initialized! \n"); return;}
	    
	    //file management
	    double param[] = {(double)L, probs.p2};
	    FILE* fp = initialize_file("data/average_noise_L%d_p2_%.2f.csv", param);
	    if(!check_pointer(fp)){printf("File could not be generated! \n"); return;}
	    double T_max = 300;
	    
	    for(int i = 0; i < N; i++){
		probs.p1 = i*delta_p;
		//allocating an array, where the mean for each p1 can be saved
		double* mean_arr = (double*) malloc(sizeof(double)* iterations);
		if(!check_pointer(mean_arr)){printf("Array could not be allocated!\n"); return;}
		
		for(int j = 0; j < iterations; j++){ //iterating over the desired number of iterations and calculating the values respectively for every iteration
		    initialize(L, grid, 0 ,generator, &probs);
		    double time = 0;
		    double delta_t = 1;
		    double infection_rate_sum = 0;
		    while(time <= T_max){
		        for(int i = 0; i < pow(L,2); i++){
		            grid_step(time, L, grid, generator, &probs);
		        }
		        infection_rate_sum += infection_rate(L, grid);
		        time += delta_t;
		    }
		    //calculating average infection rate and append it to mean array
		    double infection_rate_tempMean = infection_rate_sum / T_max;
		    mean_arr[j] = infection_rate_tempMean;
		}
		double Mean = mean(mean_arr, iterations);
		double standard = standard_deviation(mean_arr, iterations);
		fprintf(fp, "%f, %f, %f\n", probs.p1, Mean, standard);
		free(mean_arr);
		}
free(grid);
fclose(fp);
return;
	
}

    
int main(int argc, char* argv[]){
	gsl_rng* generator;
	generator = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(generator, time(NULL));
	
	
	// Task 1: Goodness of random numbers
	{
		//autocorrelation
		int len = 5000;
		double* rand_auto = (double*) malloc(sizeof(double) * len);
		if(!check_pointer(rand_auto)){printf("Memory could not be allocated! \n"); return -1;}
		for(int i = 0; i < len; i++){
		    rand_auto[i] = ran_num(generator);
		}
		auto_correlation(rand_auto, len);
		free(rand_auto);
		
		//3d plot
		len = 60000;
		double* rand_3d = (double*) malloc(sizeof(double) * len);
		if(!check_pointer(rand_3d)){printf("Memory could not be allocated! \n"); return -1;}
		for(int i = 0; i < len; i++){
		    rand_3d[i] = ran_num(generator);
		}
		goodness_of_random(rand_3d, len);
		free(rand_3d);
    	}
    	
    	
	// Task 2: creation of dataset for animation
    	{
		double p1 = 0.4;
		double p2 = 0.3;
		double p3 = 0.3;
		double p4 = 0.0;
		Probabilities probs = {p1,p2,p3,p4};
        	//For even initializations just like in the mandatory tasks -> initialize; For hotspot-spreading -> initialize_hotspots
		initFunction init_func = initialize_hotspots;
		int include_vaccination = 1; // 1: Vacc. considered 0: Vacc. not considered
		double T_max = 100;
		int L = 250;
		
		Grid(*grid)[L+2] = calloc(sizeof(Grid[L+2][L+2]),1);
		if(!check_pointer(grid)){printf("System ran out of memory!\n"); return -1;}
		
		animation(L, T_max, grid, probs, include_vaccination, generator, init_func);
		free(grid);
    	}
	
	{
		// Calculating the infection rates for 16x16 grid for multiple times and and taking the mean and standarddeviation
		double p1 = 0.0;
		double p2 = 0.6;
		double p3 = 0.3;
		double p4 = 0.0;
		Probabilities probs = {p1,p2,p3,p4};
	    
		int include_vaccination = 0; // 1: True 0: False
		int L = 16;
		int iterations = 50;
		
		average_noise(L, iterations,generator, &probs);
	    
    	}
	
	gsl_rng_free(generator);
	return 0;
	
}
