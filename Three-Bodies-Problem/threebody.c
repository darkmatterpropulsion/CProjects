/**This code simulates the motion of three bodies subjected only to their reciprocal gravitational force. The hypothesis will be that the masses satisfy M1=M2=M 
and m3 << M**/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void solve(double h,double t[],double z[][5],double T);
void bashcroftstep(double h,double t[],double z[][5],double T,double approx[12][2]);
void moultonstep(double h,double t[],double z[][5],double T,double approx[12][2]);
void rk4(double h,double t[],double z[][5]);
void g(double t,double z[][5],double vect[],int ind);
int main()
{

	//M1 initial conditions
	double y10 = 0.0;
	double vy10 = 0.0;
	double x10 = 0.1;
	double vx10 = 0.0;
	
	//M2 initial conditions
	double y20 = 0.0;
	double vy20 = 0.5;
	double x20 = 1.0;
	double vx20 = 0.0;
	
	//m3 initial conditions
	double y30 = 0.0;
	double vy30 = 0.0;
	double x30 = 0.6;
	double vx30 = 0.0;
	
	double h = 1e-6; //initial step for the rk4 method
	double t0 = 0.0;
	double T = 7; //right end of the time interval in which we search for the solution
	double z[12][5] = {};
	double t[5] = {};
	t[0] = t0;
		
	//M1
	z[0][0] = y10;
	z[1][0] = vy10;
	z[2][0] = x10;
	z[3][0] = vx10;
	
	//M2
	z[4][0] = y20;
	z[5][0] = vy20;
	z[6][0] = x20;
	z[7][0] = vx20;
	
	//m3
	z[8][0] = y30;
	z[9][0] = vy30;
	z[10][0] = x30;
	z[11][0] = vx30;
	
	rk4(h,t,z);
	solve(h,t,z,T);
}

void solve(double h,double t[],double z[][5],double T)
{
	int i;
	double approx[12][2] = {};
	FILE *out;
	out = fopen("result","w");
	fprintf(out,"%lf %lf %lf %lf %lf %lf %lf\n",t[0],z[2][0],z[0][0],z[6][0],z[4][0],z[10][0],z[8][0]);
	fprintf(out,"%lf %lf %lf %lf %lf %lf %lf\n",t[1],z[2][1],z[0][1],z[6][1],z[4][1],z[10][1],z[8][1]);
	fprintf(out,"%lf %lf %lf %lf %lf %lf %lf\n",t[2],z[2][2],z[0][2],z[6][2],z[4][2],z[10][2],z[8][2]);
	fprintf(out,"%lf %lf %lf %lf %lf %lf %lf\n",t[3],z[2][3],z[0][3],z[6][3],z[4][3],z[10][3],z[8][3]);

	t[4] = t[3] + h;
	
	while(t[4] <= T)
	{	
		bashcroftstep(h,t,z,T,approx);
		moultonstep(h,t,z,T,approx);

		fprintf(out,"%lf %lf %lf %lf %lf %lf %lf\n",t[4],z[2][4],z[0][4],z[6][4],z[4][4],z[10][4],z[8][4]);
		
		for(i = 0;i < 12;i++)
		{
			z[i][0] = z[i][1];
			z[i][1] = z[i][2];
			z[i][2] = z[i][3];
			z[i][3] = z[i][4];
		}
		
		t[0] = t[1];
		t[1] = t[2];
		t[2] = t[3];
		t[3] = t[4];
		
		t[4] = t[3] + h;
	}
	
	fclose(out);
}
void rk4(double h,double t[],double z[][5])
{
	double k[4];
	double vect[12];
	double z_temp[12][5] = {};
	int i,j,r,ind;
	for(j = 1; j < 4;j++)
	{
		ind = j-1;
		for(i = 0;i < 12;i++)
		{
			g(t[j-1],z,vect,ind);
			k[0] = vect[i];
			
			for(r = 0;r < 12;r++)
			{
				z_temp[r][ind] = z[r][ind] + 0.5*h*k[0];
			} 
			
			g(t[j-1] + 0.5*h,z_temp,vect,ind);
			k[1] = vect[i];
			for(r = 0;r < 12;r++)
			{
				z_temp[r][ind] = z[r][ind] + 0.5*h*k[1];
			} 
			
			g(t[j-1] + 0.5*h,z_temp,vect,ind);
			k[2] = vect[i];
			
			for(r = 0;r < 12;r++)
			{
				z_temp[r][ind] = z[r][ind] + h*k[2];
			} 
			
			g(t[j-1] + h,z_temp,vect,ind);
			k[3] = vect[i];
			
			z[i][j] = z[i][j-1] + h/6.0*(k[0] + 2*k[1] + 2*k[2] + k[3]);
			t[j] = t[j-1] + h;

		}	
	}
}
void bashcroftstep(double h,double t[],double z[][5],double T,double approx[][2])
{
	int i;
	double vect0[12] = {},vect1[12] = {},vect2[12] = {},vect3[12] = {};
	g(t[0],z,vect0,0);
	g(t[1],z,vect1,1);
	g(t[2],z,vect2,2);
	g(t[3],z,vect3,3);
	
	for(i = 0;i < 12;i++)
	{
		z[i][4] = z[i][3] + h/24.0*(55*vect3[i] - 59*vect2[i] + 37*vect1[i] - 9*vect0[i]);
		approx[i][0] = z[i][4];
	}
	
	h = h/2;
	for(i = 0;i < 12;i++)
	{
		approx[i][1] = z[i][3] + h/24.0*(55*vect3[i] - 59*vect2[i] + 37*vect1[i] - 9*vect0[i]);
	}
	
	
}
void moultonstep(double h,double t[],double z[][5],double T,double approx[][2])
{
	int i;
	double x0,x1,x2,f0,f1;
	double vect1[12] = {},vect2[12] = {},vect3[12] = {},vect4[12] = {},vect41[12] = {};
	g(t[1],z,vect1,1);
	g(t[2],z,vect2,2);
	g(t[3],z,vect3,3);
	g(t[4],z,vect4,4);
	
	for(i = 0;i < 12;i++)
	{
		z[i][4] = approx[i][1];
	}
	g(t[4],z,vect41,4);
	
	
	for(i = 0;i < 12;i++)
	{
		x0 = approx[i][0];
		x1 = approx[i][1];
		
		f0 = approx[i][0] - z[i][3] - h/24.0*(9*vect4[i] + 19*vect3[i] -5*vect2[i] + vect1[i]);
		f1 = approx[i][1] - z[i][3] - h/24.0*(9*vect41[i] + 19*vect3[i] -5*vect2[i] + vect1[i]);
		
		x2 = x1 - (x1-x0)/(f1-f0)*f1;
		z[i][4] = x2;
	}
	
	
}
void g(double t,double z[][5],double vect[],int ind)
{
	double y1 = z[0][ind];
	double vy1 = z[1][ind];
	double  x1 = z[2][ind];
	double vx1 = z[3][ind];
	
	double y2 = z[4][ind];
	double vy2 = z[5][ind];
	double x2 = z[6][ind];
	double vx2 = z[7][ind];
	
	double y3 = z[8][ind];
	double vy3 = z[9][ind];
	double x3 = z[10][ind];
	double vx3 = z[11][ind];
	
	//M1
	vect[0] = vy1;
	vect[1] = -(y1 - y2)/pow((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2),3.0/2.0);
	vect[2] = vx1;
	vect[3] = -(x1 - x2)/pow((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2),3.0/2.0);

	//M2
	vect[4] = vy2;
	vect[5] = -(y2-y1)/pow((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1),3.0/2.0);
	vect[6] = vx2;
	vect[7] = -(x2-x1)/pow((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1),3.0/2.0);
	
	//m3
	vect[8] = vy3;
	vect[9] = -(y3-y1)/pow((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3),3.0/2.0) - (y3-y2)/pow((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3),3.0/2.0);
	vect[10] = vx3;
	vect[11] = -(x3-x1)/pow((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3),3.0/2.0) - (x3-x2)/pow((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3),3.0/2.0);	
}
