/**This codes solves the problem of the traveling salesman. We have a set of "cities" given by [xgraph,ygraph] that are initially run in the order by which they are written, with a total distance traveled by the "salesman" of calculated as a cost function to be minimized. The minimization of the cost function is done through the method of SIMULATED ANNEALING. For 36 cities almost always I can arrive to an optimum solution with my average computer. **/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

unsigned long long int seed = 0;
double rng_glc();
void init_rng();

void init_city(double xgraph[],double ygraph[],int number_of_city);
void print_tofile(double xgraph[],double ygraph[],int exchange[],int number_of_city,FILE*path);
double calc_dist(double xgraph[],double ygraph[],int exchange[],int number_of_city);
void swap (int exchange[],int shift1,int shift2);
double simulated_annealing(double xgraph[],double ygraph[],int exchange[],int max_iteration,double *T,double T_min,double cooling_factor,int number_of_city,int *n,FILE *out);
void reset(int exchange[],int number_of_city);

int main()
{
	init_rng();
	FILE * out;
	FILE *postopt;
	double T,T_min,T_copy,best_T,cooling_factor,best_energy,energy;
	int number_of_city,max_iteration,n,best_n,run,best_run,count;
	char filename1[30] = "postopt"; //will contain the optimum path
	char filename2[30] = "cost"; // values of the cost function (total distance)
	number_of_city = 36;
	/**Test point, the traveling path is in the order by which they are presented **/
		double xgraph[36] = {0.202693,0.266395,0.824269,0.533117,0.029325,0.103711,0.719477,0.255665,0.682947,0.166787,0.720317,
				0.107313,0.200079,0.153476,0.838256,0.572539,0.368501,0.337629,0.188067,0.466100,0.354024,0.991307,0.138793,0.250223,
				0.943528,0.553530,0.654419,0.315718,0.266722,0.288295,0.257462,0.869597,0.282050,0.283195,0.558280,0.202693};
		double ygraph[36] = {0.66676,0.29748,0.49048,0.10192,0.86547,0.07554,0.25175,0.95339,0.28881,0.19389,0.37598,0.61154,
				0.72387,0.47163,0.56170,0.66028,0.39110,0.52401,0.84170,0.73495,0.07675,0.90445,0.69050,0.49740,0.86806,0.17627,0.82080,
				0.26627,0.79361,0.37283,0.15653,0.32357,0.41091,0.65368,0.01542,0.66676};
	int exchange[36]; // array in which exchange are conserved
	reset(exchange,number_of_city);//inizializzo l'array exchange
	printf("Number of run to execute:");
	scanf("%d",&run); //number of repeated run
	printf("\nMax number of iterations each run: ");
	scanf("%d",&max_iteration); //maximum number of iteration
	printf("\nInitial temperature: ");
	scanf("%lf",&T_copy); //Initial temperature
	printf("\nMin. temperature: ");
	scanf("%lf",&T_min); //Min. temperature
	printf("\nCooling factor: ");
	scanf("%lf",&cooling_factor); //the cooling factor, that is the "speed" at which the cooling has to happen.
	count = 0;

	best_energy = calc_dist(xgraph,ygraph,exchange,number_of_city); //distance calculation for the intial configuration
	while(count < run)
	{
		T = T_copy;
		printf("***********RUN %d*************\n",count);
		n=0;
		snprintf(filename1,20,"%d",count);
		strcat(filename1,"postopt");
		postopt = fopen(filename1,"w");

		snprintf(filename2,20,"%d",count);
		strcat(filename2,"cost");
		out = fopen(filename2,"w");
		//simulate annealing step with temperature T.
		energy = simulated_annealing(xgraph,ygraph,exchange,max_iteration,&T,T_min,cooling_factor,number_of_city,&n,out);
		printf("energy = %.8lf\n",energy);
		printf("Temperature rached = %.8lf in %d iteration\n",T,n);
		print_tofile(xgraph,ygraph,exchange,number_of_city,postopt);
		//store the parameters of the best run, to present them as the output of our program
		if(energy<best_energy)
		{
			best_energy = energy;
			best_T = T;
			best_n = n;
			best_run = count;
		}

		count++;
		reset(exchange,number_of_city);
		fclose(postopt);
		fclose(out);
		sleep(2); // Restart with a new temperature, by waiting 2 seconds to ensure that we generate a NEW random number (the time is the seed). 
	}
	printf("*******RUN RESULTS********\n\n");
	printf("Best run = %d\n",best_run);
	printf("Energy of the run = %lf\n",best_energy);
	printf("Temperature reach by the run = %lf\n",best_T);
	printf("Number of iteration = %d\n\n",best_n);
}

double simulated_annealing(double xgraph[],double ygraph[],int exchange[],int max_iteration,double *T,double T_min,double cooling_factor,int number_of_city,int *n,FILE *out)
{
	double d0,d1;
	int shift1,shift2;

	d0 = calc_dist(xgraph,ygraph,exchange,number_of_city);
	while((*n) < max_iteration && (*T) >= T_min)
	{
		//generation of two random number to be used as shift location
		shift1 = 1 + (abs(rng_glc()*100) % (number_of_city-2));
		shift2 = 1 + (abs(rng_glc()*100) % (number_of_city-2));
//making the shift
		swap(exchange,shift1,shift2);
//calculation of the total distance of the new configuration
		d1 = calc_dist(xgraph,ygraph,exchange,number_of_city);
//cheking if the swap should be accepted, has requested by the simulated annealing method.
		if(d1 > d0)
		{
			if(rng_glc() <= exp(-(d1-d0)/(*T)))
			{
				d0 = d1;
				*T = (*T)*cooling_factor;
			}else
			{
				swap(exchange,shift2,shift1); //reject the swap
			}
		}else
		{
			d0 = d1;
		}
		fprintf(out,"%d %.10lf %.10lf \n",*n,d0,*T);
		(*n)++;
	}
	return d0;
}

double calc_dist(double xgraph[],double ygraph[],int exchange[],int N)
{
	//The distance that we use is the euclidean distance. 
	int k;
	double distance = 0;
	for(k = 1;k < N;k++)
	{
		distance += sqrt(pow(xgraph[exchange[k]] - xgraph[exchange[k-1]],2) + pow(ygraph[exchange[k]] - ygraph[exchange[k-1]],2));
	}

	return distance;
}

void swap (int exchange[],int shift1,int shift2)
{
	double temp = 0;
	temp = exchange[shift1];
	exchange[shift1] = exchange[shift2];
	exchange[shift2] = temp;
}

void print_tofile(double xgraph[],double ygraph[],int exchange[],int number_of_city,FILE* path)
{
	int k;
	for(k = 0;k < number_of_city;k++)
	{
		fprintf(path,"%lf %lf\n",xgraph[exchange[k]],ygraph[exchange[k]]);
	}
}
void reset(int exchange[],int number_of_city)
{
	int k;
	for(k = 0;k < number_of_city;k++)
	{
		exchange[k] = k;
	}
}



//random number generator. As required by the exercise, I cannot use the random number generator of the C language.
double rng_glc()
{
	long long int t = 16807*(seed % 548781581296767) - 12838*(seed/548781581296767);

	if(t > 0)
	{
		seed = t;
	}
	else
	{
		seed = 9223372036854775807 + t;
	}
	return (double)seed/9223372036854775808;
}


void init_rng()
{
	seed = time(NULL);
}
