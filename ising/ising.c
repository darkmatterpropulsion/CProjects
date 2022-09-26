/**This program implement the method of SIMULATED ANNEALING to solve the 2-dimensional ising model.**/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

unsigned long long int seed = 0;
double rng_glc();
void init_rng();

void reset(int N,int A[][N]);
void print_matrix(int N,int A[][N]);
double calc_energy(int N,int A[][N]);
double calc_magn(int N,int A[][N]);
void find_neighbors(int N,int neighbors[],int shiftx,int shifty);
void ising(int N,int A[][N],int number_of_configuration,int PAD,double T,double *avg_energy,double *mag);
int main()
{
	init_rng();
	FILE *out;
	int N,A[30][30],number_of_configuration,PAD;
	double mag,avg_energy,T;
	N = 30; //Grid NxN
	T = 15; // Initial temperature
	number_of_configuration = 500000; // Number of configurations to try for each run.
	PAD = 1000; //Number of initial iterations, serves to reduce the correlations between various configurations. 
	out = fopen("res","w");

	while(T >= 0.1)
	{
		printf("T = %lf \n",T);
		avg_energy = 0; // average energy per site
		mag = 0; //average mangetization per site
		reset(N,A); //reset A to initial configuration

		ising(N,A,number_of_configuration,PAD,T,&avg_energy,&mag); //solving the configuration

		//averaging the quantities by keeping in mind that we sample at each 12 points of the configurations (to avoid strong correlation between configuration)
		avg_energy = avg_energy/(number_of_configuration/12);
		mag = mag/(number_of_configuration/12);

		fprintf(out,"%lf %lf %lf\n",T,avg_energy,mag);

		T = T - 0.1; //0.1 is the reduction step.
	}
	fclose(out);
}
void ising(int N,int A[][N],int number_of_configuration,int PAD,double T,double *avg_energy,double *mag)
{
	int shiftx,shifty,k,neighbors[4],de;
	double prob[2],temp;
	k = 0;

	//probability lookup table
	prob[0] = exp((-8)/T);
	prob[1] = exp((-4)/T);



	while(k < number_of_configuration+PAD)
	{
		//number between 0 and N are the spin to invert
		shiftx = abs(1000*rng_glc()) % N;
		shifty = abs(1000*rng_glc()) % N;

		//Find the coordinates of the nearest neighbors, keeping in mind the periodic boundary conditions. 
		find_neighbors(N,neighbors,shiftx,shifty);
		//changes in energy between the two configurations
		de = 2*A[shiftx][shifty]*(A[neighbors[0]][shifty] + A[neighbors[1]][shifty] + A[shiftx][neighbors[3]] + A[shiftx][neighbors[2]]);

		if(de <= 0)
		{
			A[shiftx][shifty] = -A[shiftx][shifty]; //I accept the swap

		}else
		{
			if(A[shiftx][shifty] > 0)
				temp = prob[0];
			else
				temp = prob[1];

			if(rng_glc() <= temp)
			{
				A[shiftx][shifty] = -A[shiftx][shifty]; //Conditionally accepting the swap based on the probability of it.
			}
		}
		if(k > PAD && k % 12 == 0) //each 12 sample the energy and the magnetization
		{
			(*avg_energy) += calc_energy(N,A);
			(*mag) += calc_magn(N,A);
		}
		k++;
	}

	//energy and magnetization PER SITE
	(*avg_energy) /= (N*N);
	(*mag) /= (N*N);
}
void find_neighbors(int N,int neighbors[],int shiftx,int shifty)
{
	//neighbors[0] top
	//neighbors[1] bottom
	//neighbors[2] left
	//neighbors[3] right

	neighbors[0] = shiftx-1;
	neighbors[1] = shiftx+1;
	neighbors[2] = shifty-1;
	neighbors[3] = shifty+1;

	if(shiftx == 0)
	{
		neighbors[0] = N-1;
	}else if(shiftx == N-1)
	{
		neighbors[1] = 0;
	}
	if(shifty == 0)
	{
		neighbors[2] = N-1;

	}else if(shifty == N-1)
	{
		neighbors[3] = 0;
	}


}

double calc_magn(int N,int A[][N])
{
	int i,j;
	double mag=0;;
	for(i=0;i < N;i++)
	{
		for(j=0;j < N;j++)
		{
			mag += A[i][j];
		}
	}
	return (double)mag;
}
double calc_energy(int N,int A[][N])
{
	int i,j;
	double energy = 0;
	int top,bottom,left,right;
	for(i = 0;i < N;i++)
	{
		top = i-1;
		bottom = i+1;

		if(i == 0)
			top = N-1;
		else if(i == N-1)
			bottom = 0;

		for(j = 0;j < N;j++)
		{
			left = j-1;
			right = j+1;

			if(j == 0)
				left = N-1;
			else if(j == N-1)
				right = 0;

			energy += (double)A[i][j]*(A[bottom][j] + A[i][right] + A[i][left] + A[top][j]);
		}
	}
	return -energy;
}
void reset(int N,int A[][N])
{
	int i,k;
	for(k = 0;k < N;k++)
	{
		for(i = 0;i < N;i++)
		{
			A[k][i] = 1;
		}
	}
}

void print_matrix(int N,int A[][N])
{
	int i,k;
	for(k = 0;k < 10;k++)
	{
		for(i=0;i< 10;i++)
		{
			printf("%d ",A[k][i]);
		}
		printf("\n");
	}
}
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
