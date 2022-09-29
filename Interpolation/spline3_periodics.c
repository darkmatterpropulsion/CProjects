/**This program implements the cubic splines with periodics boundary conditions. The implementation is almost identical to the natural splines, except for the definition of the coefficients.**/
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

void LU(int N,int M,double A[][M]);
void print_matrix(int N,int M,double A[][M]);
int  find_index(double xval,int N, double xdata[]);
void backsub(int N,double A[][N],double b[],double M[]);
void fwsubWithUnitaryDiag(int N,double A[][N],double y[],double delta[]);
void init_interval(double a_int,double b_int,double x[],int fpoints);
void build_matrix(int N, double A[][N],double l[],double u[]);
void build_vect(double l[],double y[],double u[],double h[],int N,double xdata[],double ydata[]);
void build_coefficient(double a[],double b[],double c[],double d[],int N,double h[],double M[],double ydata[]);

int main(void)
{
	double a_int = 2;
	double b_int = 8;
	int N = 7;
	int i;
	int fpoints = 1000;
	double xdata[7] = {2,3,4,5,6,7,8}; // nodi
	double ydata[7] = {1.4838634,1.2336425,0.81700126,0.48566875,0.28256609,0.18395886,0.16947156};
	double *l = (double*)calloc(N,sizeof(double));
	double *u = (double*)calloc(N,sizeof(double));
	double *y = (double*)calloc(N,sizeof(double));
	double *M = (double*)calloc(N,sizeof(double));
	double *h = (double*)calloc(N,sizeof(double));
	double *delta = (double*)calloc(N,sizeof(double));

	double *a = (double*)calloc(N,sizeof(double));
	double *b = (double*)calloc(N,sizeof(double));
	double *c = (double*)calloc(N,sizeof(double));
	double *d = (double*)calloc(N,sizeof(double));
	double *f = (double*)calloc(fpoints,sizeof(double));
	double A[6][6];
	double x[fpoints];


	/*intializing the discretized interval in fpoints points*/
	init_interval(a_int,b_int,x,fpoints);
	//building the l y u h vectors necessary to the costruction of the splines
	build_vect(l,y,u,h,N,xdata,ydata);

	// building the coefficient matrix A*M = y
	build_matrix(N-1,A,l,u);

	print_matrix(N-1,N-1,A);

	//Factorizing A in L, a lower triangular matrix,  and U, an upper triangular matrix. 
	//I save both matrices in A, just to saves some space. 
	LU(N-1,N-1,A);

	printf("\n\n");

	print_matrix(N-1,N-1,A);

	//making a forward substitution to solve L*Delta=y, where Delta = U*M
	fwsubWithUnitaryDiag(N-1,A,y,delta);


	//Once we found the matrix Delta, we find M with a backward substitution.
	backsub(N-1,A,delta,&M[1]);
	M[0] = M[N-1];
	for(i = 0;i < N;i++)
		printf("M[%d] = %lf\n",i,M[i]);



	//building polynomial's coefficients
	build_coefficient(a,b,c,d,N,h,M,ydata);
	
	//output of the results to file.
	FILE *out;
	out = fopen("sp","w");
	int index = 0;
	for(i=0;i < fpoints;i++)
	{
		index = find_index(x[i],N,xdata);
		f[i] = a[index] + b[index]*(x[i] - xdata[index]) + c[index]*pow(x[i] - xdata[index],2) + d[index]*pow(x[i] - xdata[index],3);
		fprintf(out,"%lf %lf\n",x[i],f[i]);
	}
	
	fclose(out);
	free(l);
	free(u);
	free(y);
	free(M);
	free(h);
	free(a);
	free(b);
	free(c);
	free(d);
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

void build_vect(double l[],double y[],double u[],double h[],int N,double xdata[],double ydata[])
{
	int j;
	for(j = 0;j <= N-2;j++)
	{
		h[j] = xdata[j+1] - xdata[j];
		printf("h[%d] = %lf \n",j,h[j]);
	}
	
	for(j = 0;j <= N-2;j++)
	{
		if(j == N-2)
		{
			u[j] = h[j]/(h[j] + h[0]);
			l[j] = h[0]/(h[j] + h[0]);
			y[j] = ((ydata[0] - ydata[j+1])/h[0] - (ydata[j+1] - ydata[j])/h[j])*6/(h[j] + h[0]);
			printf("u[%d] = %lf l[%d] = %lf  y[%d] = %lf \n",j,u[j],j,l[j],j,y[j]);
			continue;
		}
		u[j] = h[j]/(h[j] + h[j+1]);
		l[j] = h[j+1]/(h[j] + h[j+1]);
		y[j] = ((ydata[j+2] - ydata[j+1])/h[j+1] - (ydata[j+1] - ydata[j])/h[j])*6/(h[j] + h[j+1]);

		printf("u[%d] = %lf l[%d] = %lf  y[%d] = %lf \n",j,u[j],j,l[j],j,y[j]);
	}


}

void fwsubWithUnitaryDiag(int N,double A[][N],double y[],double delta[])
{
	int i,j;
	for(i = 0;i < N;i++)
	{
		delta[i] = y[i];

		for(j=0;j <= i-1;j++)
			delta[i] = delta[i] - A[i][j]*delta[j];
	}

}

void build_matrix(int N, double A[][N],double l[],double u[])
{

	printf("----Build matrix\n");
	int i,p,j,q;
	p = 0;
	q = 0;
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
				A[i][j] = 2;
				continue;
			}else if(j == i+1)
			{
				A[i][j] = u[q];
				q++;
				continue;
			}else if(i == 0 && j == N-1)
			{
				A[i][j] = l[0];
				continue;
			}else if(i == N-1 && j == 0)
			{
				A[i][j] = u[0];
				continue;
			}
			A[i][j] = 0;
		}
	}
}

void build_coefficient(double a[],double b[],double c[],double d[],int N,double h[],double M[],double ydata[])
{
	int k;
	for(k = 1; k <= N-1;k++)
	{
		a[k-1] = ydata[k-1];
		b[k-1] = (ydata[k] - ydata[k-1])/h[k-1] - (double)1/6*(M[k] + 2*M[k-1])*h[k-1];
		c[k-1] = (double)M[k-1]/2;
		d[k-1] = (double)(M[k] - M[k-1])/(6*h[k-1]);
		printf("a[%d] = %lf b[%d] = %lf c[%d] = %lf d[%d] = %lf \n",k-1,a[k-1],k-1,b[k-1],k-1,c[k-1],k-1,d[k-1]);
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

