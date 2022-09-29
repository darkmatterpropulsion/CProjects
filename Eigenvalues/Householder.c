/**This program implements the method of householder fot the tridiagonalizzation of a matrix, with subsequent use of a recursive formula for the solution of the characteristic polynomial of the matrix. In output there will be the roots of the characterisitc polynomial, that is the eigenvalues. **/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

void print_matrix(int N,double A[][N]);
void print_array(int N,double a[]);
void householder(int N,double A[][N]);
double buildsigma(int N,double A[][N],int i);
void calc_p(int N,double A[][N],double u[],double p[]);
double calc_k(int N,double u[],double p[]);
void reset(int N,double arr[]);
double find_norm(int N,double A[][N]);
int number_of_root(int N,int count[]);
double polynomial(int N,double A[][N],double x,int i,int count[]);
void find_root(int N,double A[][N],double regn[],double xs,double xd,int count[]);
void find_rootv2(int N,double A[][N],double regn[],double autoval[],double xs,double xd,int count[]);
double bisezione(int N,double A[][N],double a,double b,double toll,int count[]);


int main()
{
	int i,N = 6;
	double max;
	double A[6][6] = {{1.365,4.556,1.201,7.050,4.336,8.460},{4.556,7.343,7.797,2.840,7.662,0.896},{1.201,7.797,5.850,3.847,7.598,0.414},{7.050,2.840,3.847,5.317,2.822,8.794},{4.336,7.662,7.598,2.822,5.911,2.350},{8.460,0.896,0.414,8.794,2.350,1.708}}; //test matrix

	
	double *regn = (double *)calloc(sizeof(double),2*N);
	int *count = (int*)calloc(sizeof(int),N+1);
	double *autoval = (double *)calloc(sizeof(double),N);
	
	print_matrix(N,A);
	printf("\n\n");

	householder(N,A); //apply householder to the matrix
	print_matrix(N,A);
	
	max = find_norm(N,A); // find matrix norm to be used as boundary for the search of the zeroes. 
	/**This routine is base on a theoretical analysis that gives an upper bound to the interval for the search
	of the characteristic's polynomial roots. I implemented it in a recursive way. **/
	find_rootv2(N,A,regn,autoval,-max,max,count);
	
	
	for(i = 0;i < N;i++)
	{
		printf("rad == %lf\n",autoval[i]);
	}
	free(regn);
	free(count);
}

void find_rootv2(int N,double A[][N],double regn[],double autoval[],double xs,double xd,int count[])
{
	double x;
	int nroot_s,nroot_d,nroot_x;
	static int index = 0;
	
	polynomial(N,A,xs,N,count);
	nroot_s = number_of_root(N,count);
	
	polynomial(N,A,xd,N,count);
	nroot_d = number_of_root(N,count);
	
	x = (xs + xd)/2;
	polynomial(N,A,x,N,count);
	nroot_x = number_of_root(N,count);
	
	if(nroot_s - nroot_x == 1)
	{
		autoval[index] = bisezione(N,A,xs,x,1e-6,count);
		index++;
	}else if(nroot_s - nroot_x > 1)
	{
		find_rootv2(N,A,regn,autoval,xs,x,count);
	}
	
	if(nroot_x - nroot_d == 1)
	{
		autoval[index] = bisezione(N,A,x,xd,1e-6,count);
		index++;
	}else if(nroot_x - nroot_d > 1)
	{
		find_rootv2(N,A,regn,autoval,x,xd,count);
		
	}
	
}

void find_root(int N,double A[][N],double regn[],double xs,double xd,int count[])
{
	int index = 0,nroot_prec,nroot_new = 1;
	double tn = xs,h = 0.1; //step
	
	polynomial(N,A,tn,N,count);
	nroot_prec = number_of_root(N,count);
	
	while(tn <= xd || nroot_new != 0) //cycle until the end of the interval or until there are no remaining roots.
	{
		tn += h;
		polynomial(N,A,tn,N,count);
		nroot_new = number_of_root(N,count);
		if(nroot_new == nroot_prec-1)
		{
			regn[index] = tn-h;
			regn[index + N] = tn;
			nroot_prec = nroot_new;
			index++;
			h = 0.1; //return to the previous step
		}else if(nroot_new < nroot_prec - 1)
		{
			tn = tn-h;  //went too far
			h = h/2;   //halving the step and retry
		}
		
	}
}

int number_of_root(int N,int count[]) //this routing count the number of roots of the polynomials
{
	int i,changes = 0;
	for(i = 0;i < N;i++)
		if(count[i] == count[i+1])
			changes++;
	return changes;
}

double bisezione(int N,double A[][N],double a,double b,double toll,int count[]) //uso bisezione.
{	
	double c,fc;
	double r0 = a;
	double r1 = b;
	double fa = polynomial(N,A,a,N,count);
	double fb = polynomial(N,A,b,N,count);

	while(fabs(r0 - r1) > toll)
	{
		c = (b + a)/2;
		fc = polynomial(N,A,c,N,count);
		r0 = r1;
		if(fa*fc < 0)
		{
			b = c;
			fb = fc;
		}
		else if(fc*fb < 0)
		{
			a = c;
			fa = fc;
		}
		r1 = c;
	}

	return r1;
}
double find_norm(int N,double A[][N])
{
	int i,j;
	double max = 0,max_new;
	for(i = 0;i < N;i++)
	{
		max_new = 0;
		for(j = 0;j < N;j++)
		{
			max_new += fabs(A[i][j]);
			if(max_new > max)
			{
				max = max_new;
			}
		}
	}
	return max;	
}

void householder(int N,double A[][N])
{
	int i,j,r,c;
	double K,epsilon = 0;

	double* u = (double*)calloc(sizeof(double),N);
	double* p = (double*)calloc(sizeof(double),N);
	double* q = (double*)calloc(sizeof(double),N);
	
	for(i = N-1;i >= 2;i--)
	{
		K = 0;
		epsilon = 0;
		
		reset(N,u);
		reset(N,p);
		reset(N,q);
		
		
		for(r = 0;r < i;r++)
			epsilon += fabs(A[i][r]);
		if(epsilon == 0.0)
			continue;
		else
		{
			for(j = 0;j < i;j++)
			{
				if(j != i-1)
					u[j] = A[i][j]/epsilon;
				else
				{
					if(A[i][j] >= 0)
						u[j] = A[i][j]/epsilon + buildsigma(N,A,i)/epsilon;
					else if(A[i][j] < 0)
						u[j] = A[i][j]/epsilon - buildsigma(N,A,i)/epsilon;
				}
				
			}
			
			calc_p(N,A,u,p);
			K = calc_k(N,u,p);
			
			for(r = 0;r < N;r++)
				q[r] = p[r] - K*u[r];
				
						
			for(r = 0;r < N;r++)
			{
				for(c = 0; c < N;c++)
				{
					A[r][c] = A[r][c] - q[r]*u[c] - u[r]*q[c];				
				}
			}

			printf("\n\n");
		}
	}
	
	free(u);
	free(q);
	free(p);
}


double calc_k(int N,double u[],double p[])
{
	int j;
	double K = 0;
	double u2;
	for(j = 0;j < N;j++)
		u2 += u[j]*u[j];
	
	for(j = 0; j < N;j++)
		K += u[j]*p[j];

	K = K/u2;

	return K;
}

void calc_p(int N,double A[][N],double u[],double p[])
{
	int i,j;
	double u2 = 0;
	for(j = 0;j < N;j++)
		u2 += u[j]*u[j];

	for(i = 0;i < N;i++)
	{
		for(j = 0;j < N;j++)
		{
			p[i] += A[i][j]*u[j]/(0.5*u2);
		}
	}
}

double buildsigma(int N,double A[][N],int i)
{
	int q;
	double sigma = 0;
	for(q = 0;q < i;q++)
	{
		sigma += A[i][q]*A[i][q];
	}

	sigma = sqrt(sigma);
	
	return sigma;
}
double polynomial(int N,double A[][N],double x,int i,int count[])
{
	double val = 0;
	if(i == 0)
	{
		count[i] = 1;
		return 1;
	}
	else if(i == 1)
	{
		val = A[0][0] - x;
		count[i] = val/fabs(val);
		return val;
	}
	else
	{
		val = (A[i-1][i-1] - x)*polynomial(N,A,x,i-1,count) - A[i-1][i-2]*A[i-1][i-2]*polynomial(N,A,x,i-2,count);
		count[i] = val/fabs(val);
		return val;
	}
		 
}


void reset(int N,double arr[])
{
	int i;
	for(i = 0;i < N;i++)
	{
		arr[i] = 0;
	}
}

void print_matrix(int N,double A[][N])
{	
	int i,j;
	for(i = 0;i < N;i++)
	{
		for(j = 0;j < N;j++)
		{
			printf("%.3lf ",A[i][j]);
		}
		printf("\n");
	}	
}
void print_array(int N,double a[])
{
	int i;
	for(i = 0;i < N;i++)
		printf("%.3lf\n",a[i]);
}
