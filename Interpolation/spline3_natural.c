/**This program implements the cubic splines interpolation with natural boundary conditions. **/
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

void set_up_M(int N,double M[]);
void LU(int N,int M,double A[][M]);
void print_matrix(int N,int M,double A[][M]);
int find_index(double xval,int N, double xdata[]);
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
	double A[5][5];
	double x[fpoints];


	/*initializing the discretized interval*/
	init_interval(a_int,b_int,x,fpoints);
	//bulding the vector l y u h necessary for the spline costruction
	build_vect(l,y,u,h,N,xdata,ydata);

	// Building matrix coefficient A*M = y
	build_matrix(N-2,A,l,u);

	print_matrix(N-2,N-2,A);

	//Factorizing A in L, a lower triangular matrix,  and U, an upper triangular matrix. 
	//I save both matrices in A, just to saves some space. 
	LU(N-2,N-2,A);

	printf("\n\n");

	print_matrix(N-2,N-2,A);

	//forward surbstitution to solve L*Delta=y, where Delta = U*M
	fwsubWithUnitaryDiag(N-2,A,y,delta);


	//Once Delta matrix has been found, I find matrix M with a backward substitution.
	backsub(N-2,A,delta,M);

	set_up_M(N,M);

	for(i = 0;i < N;i++)
	{
		printf("M[%d] = %lf\n",i,M[i]);
	}


	//build the polynomials coefficient, based on the theoretical formulae
	build_coefficient(a,b,c,d,N,h,M,ydata);
	

	//write the found spline to file.
	FILE *out;
	out = fopen("sp","w");
	int index = 0;
	for(i=0;i < fpoints;i++)
	{
		index = find_index(x[i],N,xdata);
		f[i] = a[index] + b[index]*(x[i] - xdata[index]) + c[index]*(x[i] - xdata[index])*(x[i] - xdata[index]) + d[index]*(x[i] - xdata[index])*(x[i] - xdata[index])*(x[i] - xdata[index]);
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
	
	for(j = 0;j <= N-3;j++)
	{
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
