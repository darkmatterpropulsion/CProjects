/** This program implement the polinomial lagrange interpolation**/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void init_interval(int fpoints, double x[],double a,double b);
void write_res(double* f,int fpoints, double *x,char filename[],double a,double b);
void chebyshev(double nodes[],double val[],int Nnodes);
void equi_spaced(double nodes[],double val[], int Nnodes,double a,double b);
void lagrange_interpolation(double nodes[],double val[], double x[], int Nnodes,int fpoints, double f[]);
int find_index(int fpoints,double x[],double val,double a,double b);
void set_to_number(double *arr,int len,int number);
double find_max(double f[],double x[],int fpoints,double ln);
void esercizio_1();
void esercizio_2();
void esercizio_3();

int main(void)
{
	//esercizio_1();
	//esercizio_2();
	//esercizio_2();
}

void esercizio_1()
{
/**This exercise requires to calculate the difference between the true function and the interpolating polynomial. With the calculated difference we should verify
the rest theorem.  **/
	int k;
	int fpoints = 1000;
	double x[fpoints]; //contains the points of the discretized interval.
	double a = 2;
	double b = 8;
	double nodes[4] = {2,4,6,8};
	double val[4] = {1.4838634,0.81700126,0.28256609,0.16947156}; //points of the interpolated functions at nodes[]
    double * f = (double*)malloc(sizeof(double)*(fpoints+1));
    double * reminder = (double*)malloc(sizeof(double)*(fpoints+1));

    /*discretization of the a,b interval and storing the point in the array x[] */
    init_interval(fpoints,x,a,b);

    /*lagrange interpolation on the nodes stored in the array "nodes"
     with values stroged in "val". The function need also the number of nodes
     that I want to use */

     lagrange_interpolation(nodes,val,x,4,fpoints,f);

    /*write data to a file to be used in Gnuplot.*/

    write_res(f,fpoints,x,"data",a,b);

    //This calculates the remainder between the true function f[] and the polinomial interpolation
    for(k=0;k <= fpoints;k++)
    {
    	reminder[k] = f[k] - (x[k]*x[k]*exp(1-x[k]) + 1/((11-x[k])*(11-x[k])));
    }
    //Writing the result in the file "reminder"
    write_res(reminder,fpoints,x,"reminder",a,b);

    /*Calculating the max of the remainder in the interpolation interval*/
    double max = 0;
    int ind = 0;
    for(k = 0;k <= fpoints; k++)
    {
    	if(fabs(reminder[k]) > max)
    	{
    		max = fabs(reminder[k]);
    		ind = k;
    	}

    }
    printf("reminder max at x[%d] = %lf with value %lf\n",ind,x[ind],max);

    free(f);
    free(reminder);
}

void esercizio_2()
{
/**This exercise requires to verify that the error from the interpolation goes as the square of the interpolating interval h (if we use 2 nodes) or the cube (if we use 3 nodes).**/
	int i,k,ind;
	double a = -1;
	double b = 1;
	int last_node = 0;
	int fpoints = 100000;
	double x[fpoints];
	double nodes[100];
	double val[100];
	char filename[20] = "data";
	double * f = (double*)malloc(sizeof(double)*(fpoints+1));
	double max = 0;
	double h = 0;
	double N = 0;

	//Interval discretization
	init_interval(fpoints,x,a,b);

	//Cycling on 10 to 100 nodes, to see how the error changes.
	for(i = 10;i <= 100; i +=10)
	{
		//initializing to 0 the array that contains the values of the interpolating polinomial
		set_to_number(f,fpoints + 1,0);
		//generating equispatial nodes in the interval [a,b[
		equi_spaced(nodes,val,i,a,b);
		//making an interpolation each 2 points
		for(k = 0;k < (i-1);k++)
		{
				ind = find_index(fpoints,x,nodes[k],a,b); //find the index of the value in x[] that is closer to nodes[k]
				lagrange_interpolation(&nodes[k],&val[k],&x[ind],2,fpoints/i,&f[ind]);  /**interpolate by lagrange, NOTICE: the same function lagrande_interpolation 
													is used for both two and three nodes interpolation. I had to find a clever way to do it.**/
				max += find_max(&f[ind],&x[ind],fpoints/i,nodes[k+1]); /*calculating the sum of the maxima in each interval and calculating their mean*/
		}
		/*writing the results in a file named isquared where i is the number of nodes
		 * the execution of this exercise generates many files*/
		snprintf(filename,10,"%d",i);
		strcat(filename,"squared");
		/**Find the last node intex so that write_res knows when to stop **/
		last_node = find_index(fpoints,x,nodes[i-1],a,b);
		write_res(f,fpoints,x,filename,a,x[last_node]);
		//find the step h
		h = (double)(nodes[1] - nodes[0]);
		//find the number of interval for the interpolation
		N = (double)(nodes[i-1] - nodes[0])/h;
		max = max/N; //making the mean of the preceeding points.
		printf("%lf %lf\n",h,max/(h*h));

	}

	printf("\n\n");

	//same as before, this times with groups of three nodes.
	for(i = 11;i <= 95; i +=12)
	{
		set_to_number(f,fpoints + 1,0);
		equi_spaced(nodes,val,i,a,b);
		for(k = 0;k < (i-2);k+=2)
		{
				ind = find_index(fpoints,x,nodes[k],a,b);
				lagrange_interpolation(&nodes[k],&val[k],&x[ind],3,2*fpoints/i,&f[ind]);
				max += find_max(&f[ind],&x[ind],2*fpoints/i,nodes[k+2]);
		}
		snprintf(filename,10,"%d",i);
		strcat(filename,"cubic");
		last_node = find_index(fpoints,x,nodes[i-1],a,b);
		write_res(f,fpoints,x,filename,a,x[last_node]);

		h = (double)(nodes[2] - nodes[0]);
		N = (nodes[i-1] - nodes[0])/h;
		max = max/N;
		printf("%lf %lf\n",h,max/(h*h*h));
	}
}

void esercizio_3()
{
	int i;
	int last_node = 0;
	double a = -1;
	double b = 1;
	int fpoints = 1000;
	double x[fpoints];
	double nodes[100];
	double val[100];
	char filename[20] = "data";
	double * f = (double*)malloc(sizeof(double)*(fpoints+1));
	//initializing the interval
	init_interval(fpoints,x,a,b);

	//iteration between 4 and 20 nodes at step of 2
	for(i = 4;i <= 20; i += 2)
	{
		/*setting to 0 the arrays that will contain the values of the interpolating 
		polynomial at the points of the discretized interval. */
		set_to_number(f,fpoints+1,0);
		/*generating equispatial nodes between a=-1 and b =1*/
		equi_spaced(nodes,val,i,a,b);
		/*lagrande interpolation*/
		lagrange_interpolation(nodes,val,x,i,fpoints,f);
		/*writing the results in tabular form
		 * the files have the denomination iequi where
		 * i is the number of nodes*/
		snprintf(filename,20,"%d",i);
		strcat(filename,"equi");
		/**Find the last node intex so that write_res knows when to stop **/
		last_node = find_index(fpoints,x,nodes[i-1],a,b);
		write_res(f,fpoints,x,filename,a,x[last_node]);
	}

	//come con il ciclo precedente, sola che ora genero i nodi secondo chebyshev
	for(i = 4;i <= 20; i += 2)
		{
			set_to_number(f,fpoints+1,0);
			chebyshev(nodes,val,i);
			lagrange_interpolation(nodes,val,x,i,fpoints,f);
			snprintf(filename,20,"%d",i);
			strcat(filename,"cheby");
			write_res(f,fpoints,x,filename,a,b);
		}


}
void lagrange_interpolation(double nodes[], double val[], double x[], int Nnodes, int fpoints,double f[])
{

    double *lambda = (double *)malloc(sizeof(double)*(Nnodes+1));
    double * theta = (double *)calloc(Nnodes,sizeof(double));
    double phi;
    int i,j,k,r;
    int div_flag = 0;

    set_to_number(theta,Nnodes,1);

    for(i=0;i < Nnodes;i++)
    {
    	for(j=0;j < Nnodes;j++)
    	{
    		if(j == i)
    			continue;
    		else
    			theta[i] *= nodes[i] - nodes[j];


    	}
    }



    for (k = 0; k <= fpoints; k++)
    {
    	div_flag = 0;
    	if(fabs(f[k]) != 0 && k == 0)
    		continue;

    	set_to_number(lambda,Nnodes,1);

    	for(r = 0; r < Nnodes;r++)
    	{
    		if(fabs(x[k] - nodes[r]) < 0.0000001)
    			div_flag = 1;
    	}

    	if(div_flag == 0)
    	{
    		phi = 1;
			for(i = 0; i < Nnodes;i++)
			{
				f[k] += val[i]/(theta[i]*(x[k] - nodes[i]));
				phi *= x[k] - nodes[i];
			}
				f[k] *= phi;
    	}else
    	{
    		for(i = 0;i < Nnodes;i++)
    		{
    			phi = 1;
    			for(j = 0;j < Nnodes;j++)
    			{
    				if(i==j)
    					continue;
    				else
    					phi *= x[k] - nodes[j];
    			}
    			lambda[i] = phi/theta[i];
    			f[k] += lambda[i]*val[i];
    		}
    	}
    }

free(lambda);
free(theta);

}

void init_interval(int fpoints, double x[],double a,double b)
{
    int i;
    for(i = 0; i <= fpoints; i++)
            x[i] = a + (double)i*(b-a)/fpoints;
}

void write_res(double* f,int fpoints, double *x,char filename[],double a,double b)
{
	FILE *out;
	out = fopen(filename, "w");
    int k;
    for(k = 0;k < fpoints;k++)
    {
    	if(x[k] >= a && x[k] < b)
    	{
    		fprintf(out,"%lf %lf\n",x[k],f[k]);
    	}
    }
    fclose(out);
}

void chebyshev(double nodes[],double val[],int Nnodes)
{
	int i;
	for(i=1;i <= Nnodes;i++)
	{
		nodes[i] = (double)cos((double)(i - 0.5)*acos(-1.0)/(Nnodes));
		val[i] = (double)1/(1 + 25*nodes[i]*nodes[i]);

	}
}

void equi_spaced(double nodes[],double val[], int Nnodes,double a,double b)
{
	int i;
	for(i = 0;i < Nnodes;i++)
	{
		nodes[i] = a + (double)(b-a)*i/Nnodes;
		val[i] = (double)1/(1 + 25*nodes[i]*nodes[i]);
	}

}
int find_index(int fpoints,double x[],double val,double a,double b)
{
	int i;
	int ind = 0;
	for(i = 0;i <= fpoints;i++)
		{
			if(fabs(val - x[i]) <= (double)(b-a)/fpoints)
				ind = i;
		}
	return ind;
}

void set_to_number(double *arr,int len,int number)
{
	int i;
	for(i=0;i < len;i++)
	{
		arr[i] = number;
	}
}
double find_max(double f[],double x[],int fpoints,double ln)
{
	int j;
	double u,max;
	max = 0;
	for(j = 0; j <= fpoints;j++)
	{
		if(x[j] < ln - 0.01)
		{
			u = f[j] - (double)1/(1+25*x[j]*x[j]);
			if(fabs(u) > max)
				max = fabs(u);

		}else
			break;
	}

	return max;
}

