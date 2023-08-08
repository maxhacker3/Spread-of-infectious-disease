#include "numerics_disease.h"

typedef enum{//enum contains all possible health states
    SUSCEPTIBLE = 1,
    INFECTED = 2,
    RECOVERED = 3,
    VACCINATED = 4
} Status;

typedef struct{//The datatype of every grid cell is "Grid"; Each cell assumes an element from the enum
    Status status; 
}Grid;

typedef struct{
    double p1;
    double p2;
    double p3;
    double p4;
}Probabilities;


// Printing the grid in the console. Could be used to check whether the initialization has been completed correctly
void print_grid(int L, Grid grid[][L+2]){
    for(int i = 0; i < L+2; i++){
        for(int j = 0; j < L+2; j++){
            printf("%d ", grid[i][j].status);
        }
        printf("\n");
    }
    
}

// initializes the grid with SUSCEPTIBLE, INFECTED and RECOVERED; Allocation with equal probability; 
// if include_vaccination == 1,vaccination is included in the initialization
void initialize(int L, Grid grid[][L+2], int include_vaccination, gsl_rng* generator, void* params){
    Probabilities* probs = (Probabilities*)params;
    for(int i = 1; i < L+1; i++){
        for(int j = 1; j < L+1; j++){
            //draw random number if vacc. is included
            double vacc_prob = (include_vaccination) ? ran_num(generator) : 2.0;
            //checks, if vaccination applies; else decides on random state 1,2 or 3
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

//checks the neighbours according to von-Neumann neighbourhood; every adjacent infection is counted and returned
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

// probability of infection depends on number of infected neighbours
double calc_probability(double p, int n){
    return 1 - pow(1-p,n);//binomial distribution; For explanation see the report
}

// updates the grid randomly; Core of the algorithm
void grid_step(double time, int L, Grid grid[][L+2], gsl_rng* generator, void* params){
    	Probabilities* probs = (Probabilities*)params;
    	double p1 = probs->p1; double p2 = probs->p2; double p3 = probs->p3; double p4 = probs->p4;
    
    	// Draws of a random number between 0 and L^2-1 and converting it to the corresponding grid coordinates.
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
                } else{;}break;
            
            case 2://Cell infected; p2 to recover; cell is updated if applicable
                if(rand_prob <= p2){
                    grid[i][j].status = RECOVERED;
                }else{;}break;
        
            case 3://cell recovered; p3 to become susceptible; cell is updated if applicable
                if(rand_prob <= p3){
                    grid[i][j].status = SUSCEPTIBLE;
                }else{;}break;
            
            case 4://cell is vaccinated and remains vaccinated
                grid[i][j].status = VACCINATED;
            default: break;
    }
}

//given a grid, the infectious people out of all people are determined and returned
double infection_rate(int M, Grid grid[][M+2]){
    int sum = 0;
    for(int i = 1; i < M+1;i++){
        for(int j = 0; j < M+1; j++){
            sum += (grid[i][j].status == INFECTED)? 1 : 0;
        }
    }
    return (double)(sum / pow(M,2));
}


//this function executes the main tasks; It generates the data which is used for analysing the connection between average infection rate and model-probabilities (p1 and p4 is evaluated)
void task(int L, double T_max, Grid grid[][L+2], Probabilities probability, int include_vaccination, gsl_rng* generator){
    Probabilities probs = probability;
    double delta_p = 0.02; //step-width for p1/p4
    int N = 1 / 0.02; //Total number of steps
    
    //initializing a custom file
    double param[] = {(double)L, probs.p2};
    FILE* fp_rate = initialize_file("data/infection_rate_L%d_p2_%.2f.csv", param);
    if(!check_pointer(fp_rate)){printf("The file could not be generated\n"); return;} //check correct initialization
    
    
    for(int i = 0; i < N; i++){ 
    	//case distinction
        if(include_vaccination == 0){
            probs.p1 = i*delta_p; //incrementing p1 for every iteration
        }
        else{
            probs.p4 = i*delta_p; //incrementing p4 for every iteration
        }
        initialize(L, grid, include_vaccination ,generator, &probs);
     
        double time = 0; //initial time
        double delta_t = 1; //time step-width
        double infection_rate_sum = 0; //initial infection rate
        while(time <= T_max){
            for(int i = 0; i < pow(L,2); i++){
                grid_step(time, L, grid, generator, &probs); //updating the grid randomly L^2 times
            }
            infection_rate_sum += infection_rate(L, grid);
            time += delta_t;
        }
        double infection_rate_tempMean = infection_rate_sum / T_max; //calculating the average infection rate over time
        //printing average infection rates to file
        if(include_vaccination == 0){
            fprintf(fp_rate, "%f, %f\n", probs.p1, infection_rate_tempMean);
        }
        else{
            fprintf(fp_rate, "%f, %f\n", probs.p4, infection_rate_tempMean);
        }  
    }
    fclose(fp_rate);
    return;
}




int main(int argc, char* argv[]){
	//initializing the MT19937 generator
   	gsl_rng* generator;
	generator = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(generator, time(NULL));
    	
    	
   	//System parameter
   	double T_max = 300; //grid width or height respectively
   	int L_array[] = {16,32,64,150}; //considered grid sizes
    
    
   	printf("Starting...\n");
        {
        	
            // Task 2a
            //initializing the probabilities
            double p1 = 0.0;
            double p2 = 0.3;
            double p3 = 0.3;
            double p4 = 0.0;
            Probabilities probs = {p1,p2,p3,p4}; //feed the "Probabilities" structure
            
            int include_vaccination = 0; // No vaccination
            
            for(int i = 0; i < 4; i++){
                int L = L_array[i];
                Grid(*grid)[L+2] = calloc(sizeof(Grid[L+2][L+2]),1); //allocating the grid
                if(!check_pointer(grid)){printf("The grid could not be initialized! \n"); return -1;} //check correct initialization
                
                task(L, T_max, grid, probs, include_vaccination, generator);
                printf("L = %d completed...", L);
                free(grid);
            }
            printf("Task 1/3 completed\n");
            
        }
        {
            //task 2b; setting p2 from 0.3 to 0.6; Analogous strategy
            double p1 = 0.0;
            double p2 = 0.6;
            double p3 = 0.3;
            double p4 = 0.0;
            Probabilities probs = {p1,p2,p3,p4};
            
            int include_vaccination = 0; //no vaccination
            
            for(int i = 0; i < 4; i++){
                int L = L_array[i];
                Grid(*grid)[L+2] = calloc(sizeof(Grid[L+2][L+2]),1);
                if(!check_pointer(grid)){printf("The grid could not be initialized! \n"); return -1;}
                
                task(L, T_max, grid, probs, include_vaccination, generator);
                printf("L = %d completed...", L);
                free(grid);
            }
            printf("Task 2/3 completed\n");
        }
        
    
        {
            //Task 3
            //initializing the probabilities
            double p1 = 0.5;
            double p2 = 0.5;
            double p3 = 0.5;
            double p4 = 0.0;
            Probabilities probs = {p1,p2,p3,p4};
            
            int include_vaccination = 1; // Include Vaccination
            
            for(int i = 0; i < 4; i++){
                int L = L_array[i];
                Grid(*grid)[L+2] = calloc(sizeof(Grid[L+2][L+2]),1);
                if(!check_pointer(grid)){printf("The grid could not be initialized! \n"); return -1;}
                
                task(L, T_max, grid, probs, include_vaccination, generator);
                printf("L = %d completed...", L);
                free(grid);
            }
            printf("Task 3/3 completed\n");
        }
        
      
       
    
    gsl_rng_free(generator);
    return 0;
}
