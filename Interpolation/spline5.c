/**This program implements the fifth-degree spline. **/
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

void set_up_M(int N,double M[]);
void LU(int N,int M,double A[][M]);
void print_matrix(int N,int M,double A[][M]);
int find_index(double xval,int N, double xdata[]);
void backsub(int N,double A[][N],double b[],double M[]);
void fwsubs(int N,double A[][N],double y[],double coeff[]);
void init_interval(double a_int,double b_int,double x[],int fpoints);
void build_matrix(int N, double A[][N],double l[],double u[],double t[]);
void build_vect(double l[],double y[],double u[],double h[],double t[],int N,double xdata[],double ydata[]);
void build_coefficient(double alpha[],double beta[],double gamma[],double eta[],double delta[],double sigma[],int N,double h[],double M[],double ydata[]);

int main(void)
{
	double a = 0.;
	double b = 8.;
	const int N = 40;
	int i;
	int fpoints = 1000;
	double xdata[40];
	double ydata[40];
	double *t = (double*)calloc(N,sizeof(double));
	double *l = (double*)calloc(N,sizeof(double));
	double *u = (double*)calloc(N,sizeof(double));
	double *y = (double*)calloc(N,sizeof(double));
	double *M = (double*)calloc(N,sizeof(double));
	double *h = (double*)calloc(N,sizeof(double));
	double *coeff = (double*)calloc(N,sizeof(double));

	double *alpha = (double*)calloc(N,sizeof(double));
	double *beta = (double*)calloc(N,sizeof(double));
	double *gamma = (double*)calloc(N,sizeof(double));
	double *eta = (double*)calloc(N,sizeof(double));
	double *delta = (double*)calloc(N,sizeof(double));
	double *sigma = (double*)calloc(N,sizeof(double));
	double *f = (double*)calloc(fpoints,sizeof(double));
	double A[38][38];
	double x[fpoints];

	//genereting the nodes
	for(i=0;i < 40;i++)
	{
		xdata[i] = a + (b-a)*i/N;
		ydata[i] = xdata[i]*xdata[i]*exp(1-xdata[i]) + 1/((11-xdata[i])*(11-xdata[i]));
	}

	/*initilizing the discretized interval in fpoints points.*/
	init_interval(a,b,x,fpoints);

	//building the vectors l u y h necessary to the costruction of the splines
	build_vect(l,y,u,h,t,N,xdata,ydata);

	// building the coefficient matrix A*M = y
	build_matrix(N-2,A,l,u,t);

	//Factorizing A in L, a lower triangular matrix,  and U, an upper triangular matrix. 
	//I save both matrices in A, just to saves some space. 
	LU(N-2,N-2,A);


	//making a forward substitution to solve L*COEFF = y, where COEFF = U*M
	fwsubs(N-2,A,y,coeff);


	//once the matrix COEFF has been found, i find the matrix M with a backward substitution.
	backsub(N-2,A,coeff,M);

	/*I substitute the vector of coefficients M with shifted elements, so to have M[0]=M[N-1]=0 having imposed
	  the condition of natural splines. */
	set_up_M(N,M);


	//building the coefficients of the polinomial..
	build_coefficient(alpha,beta,gamma,eta,delta,sigma,N,h,M,ydata);

/*writing the results in the file sp, così da poter visualizzare l'interpolazione con gnuplot.*/

	FILE *out;
	out = fopen("sp","w");
	int index = 0;
	for(i=0;i < fpoints;i++)
	{
		if(x[i] < xdata[39])
		{
			index = find_index(x[i],N,xdata);
			f[i] = alpha[index]*pow(x[i] - xdata[index],5) + beta[index]*pow(x[i]-xdata[index+1],5) + gamma[index]*pow(x[i]-xdata[index],3) + eta[index]*pow(x[i] - xdata[index+1],3) + delta[index]*(x[i]-xdata[index]) + sigma[index]*(x[i] - xdata[index+1]);
			fprintf(out,"%lf %lf\n",x[i],f[i]);
		}else
			break;
	}

/*close the files and free the memory.*/
	fclose(out);
	free(t);
	free(l);
	free(u);
	free(y);
	free(M);
	free(h);
	free(alpha);
	free(beta);
	free(gamma);
	free(eta);
	free(delta);
	free(sigma);
}

void backsub(int N,double A[][N],double b[],double M[])
{

	int i,j;
	for(i = 0;i < N;i++)
		M[i] = b[i];

	for(i = N-1;i >= 0; i--)
	{
		M[i] = b[i];
		for (j=i+1;j <= N-1;j++)
		{
			M[i] = M[i] - A[i][j]*M[j];
		}
		M[i] = M[i]/A[i][i];
	}

}

void print_matrix(int N,int M,double A[][M])
{
	int i,j;
	for(i = 0;i < N;i++)
	{
		for(j = 0;j < M;j++)
		{
			printf("%lf ",A[i][j]);
		}
		printf("\n");
	}
}

void LU(int N,int M,double A[][M])
{
	int i,j;
	int k =0;
	double m;
	while(k < N)
	{
		for(i=k+1;i < N;i++)
		{
			m = A[i][k]/A[k][k];
			for(j=k+1;j < N;j++)
			{
				A[i][j] = A[i][j] - m*A[k][j];
			}
			A[i][k] = m;
		}
		k++;
	}
}

void build_vect(double l[],double y[],double u[],double h[],double t[],int N,double xdata[],double ydata[])
{
	int j;
	int i;
	for(j = 0;j <= N-2;j++)
		h[j] = xdata[j+1] - xdata[j];



	for(i = 1;i <= N-2;i++)
	{
		u[i-1] = pow(h[i],3)/120 - pow(h[i],2)*(2*pow(h[i-1],2) + 3*h[i-1]*h[i] + 2*pow(h[i],2))/(36*(h[i-1] + h[i]));
		t[i-1] = pow(h[i-1],3)/30 -2*pow(h[i-1],2)*(pow(h[i-1],2) + 3*h[i-1]*h[i] + 2*pow(h[i],2))/(36*(h[i-1] + h[i])) + pow(h[i],3)/30 - 2*pow(h[i],2)*(2*pow(h[i-1],2) + 3*h[i-1]*h[i] + 2*pow(h[i],2))/(36*(h[i-1] + h[i]));
		l[i-1] = pow(h[i-1],3)/120 - pow(h[i-1],2)*(pow(h[i-1],2) + 3*h[i-1]*h[i] + 2*pow(h[i],2))/(36*(h[i-1]+ h[i]));
		y[i-1] = (ydata[i+1] - ydata[i])/h[i] - (ydata[i] - ydata[i-1])/h[i-1];

		//printf("u[%d] = %lf l[%d] = %lf t[%d] = %lf  y[%d] = %lf \n",i-1,u[i-1],i-1,l[i-1],i-1,t[i-1],i-1,y[i-1]);
	}

}

void fwsubs(int N,double A[][N],double y[],double coeff[])
{
	int i,j;
	for(i = 0;i < N;i++)
	{
		coeff[i] = y[i];

		for(j=0;j <= i-1;j++)
			coeff[i] = coeff[i] - A[i][j]*coeff[j];
	}
}

void build_matrix(int N, double A[][N],double l[],double u[],double t[])
{

	//printf("----Build matrix\n");
	int i,p,j,q,g;
	p = 0;
	q = 0;
	g = 0;
	for(i = 0;i < N;i++)
	{
		for(j = 0;j < N;j++)
		{
			if(j == i-1)
			{
				A[i][j] = l[p];
				p++;
				continue;
			}else if(j == i)
			{
				A[i][j] = t[g];
				g++;
				continue;
			}else if(j == i+1)
			{
				A[i][j] = u[q];
				q++;
				continue;
			}

			A[i][j] = 0;
		}
	}
}

void build_coefficient(double alpha[],double beta[],double gamma[],double eta[],double delta[],double sigma[],int N,double h[],double M[],double ydata[])
{
	int i;
	for(i = 0; i <= N-2;i++)
	{

		alpha[i] = M[i+1]/(120*h[i]);
		beta[i] = -M[i]/(120*h[i]);
		gamma[i] = -M[i+1]*(pow(h[i],2)+3*h[i]*h[i+1] + 2*pow(h[i+1],2))/(36*(h[i] + h[i+1]));
		eta[i] = M[i]*(2*pow(h[i],2) + 3*h[i]*h[i+1] + 2*pow(h[i],2))/(36*(h[i] + h[i+1]));
		delta[i] = ydata[i+1]/h[i] - pow(h[i],2)*gamma[i] - pow(h[i],4)*alpha[i];
		sigma[i] = -ydata[i]/h[i] - pow(h[i],2)*eta[i] - pow(h[i],4)*beta[i];
		//printf("alpha[%d] = %lf beta[%d] = %lf gamma[%d] = %lf eta[%d] = %lf  \n",i,alpha[i],i,beta[i],i,gamma[i],i,eta[i]);
	}
}

int find_index(double xval,int N, double xdata[])
{
	int i,index;
	for(i = 0;i < N;i++)
	{
		if(xval >= xdata[i] && xval < xdata[i+1])
			index = i;
	}

	return index;
}

void init_interval(double a_int,double b_int,double x[],int fpoints)
{
	int i;
	for(i = 0; i < fpoints;i++)
	{
		x[i] = a_int +  (double)(b_int-a_int)*i/fpoints;
	}

}
void set_up_M(int N,double M[])
{
	double *x = (double*)calloc(N,sizeof(double));
	int i;
	for(i=0;i <= N-2;i++)
	{
		x[i] = M[i];
	}
	for(i=0;i <= N-2;i++)
	{
		if(i == 0)
			M[0] = 0;
		else if(i == N-1)
			M[N-1] = 0;
		else
			M[i] = x[i-1];
	}

	free(x);
}
