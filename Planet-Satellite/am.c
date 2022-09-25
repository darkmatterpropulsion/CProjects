/** This code solves the motion of two bodies in space subjected to the reciprocal gravitational force. The hypothesis that M >> m imply that the gravitational pull of m on M is negligible in the forces equations. Units of GM = 1 have been used to simplify the implementations of the equations into the code. **/
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
double solve(double h,double t[],double z[][5],double T,char filename[40],double stats[2]);
void bashcroftstep(double h,double t[],double z[][5],double T,double approx[][2]);
void moultonstep(double h,double t[],double z[][5],double T,double approx[][2]);
void rk4(double h,double t[],double z[][5]);
void g(double t,double z[][5],double vect[],int ind);

int main()
{
	//Initial conditions, that gives the position and velocity of m
	double y0 = 0.0;
	double vy0 = 1.63;
	double x0 = 0.5;
	double vx0 = 0;
	
	double h = 1e-3; //time step
	double t0 = 0.0; //initial time for the simulation
	double T = 10; //final time for the simulation
	double stats[2] = {};
	char filename [40] = "result";
	double z[4][5] = {};
	double t[5] = {};
	
	z[0][0] = y0;
	z[1][0] = vy0;
	z[2][0] = x0;
	z[3][0] = vx0;
	t[0] = t0;
	
	rk4(h,t,z); //Using rk4 to find the first point to kick-start the actual used method. 
	solve(h,t,z,T,filename,stats); // Predicotr corrector based on Adams-bashcroft - Adams-multon methods.
	printf("Period == %lf Eccentricity %lf\n",stats[1],stats[0]); //period and eccentricity of the orbit of the mass m.
}

double solve(double h,double t[],double z[][5],double T,char filename[40],double stats[2])
{
	int i,avg = 0,flag = 1;
	double approx[4][2] = {};
	double periodo = 0;
	double da,dp,r = 0;
	double x0 = z[2][0];
	double y0 = z[0][0];
	
	double raverage = 0;
	double rerror = 0;
	
	FILE *out;
	out = fopen(filename,"w");
	
	fprintf(out,"%lf %lf %lf\n",t[0],z[2][0],z[0][0]);
	fprintf(out,"%lf %lf %lf\n",t[1],z[2][1],z[0][1]);
	fprintf(out,"%lf %lf %lf\n",t[2],z[2][2],z[0][2]);
	fprintf(out,"%lf %lf %lf\n",t[3],z[2][3],z[0][3]);
	da = pow(z[2][3]*z[2][3] + z[0][3]*z[0][3],1.0/2.0);
	dp = da;
	
	periodo = t[0];
	
	t[4] = t[3] + h;

	while(t[4] <= T)
	{	
		avg++;
		bashcroftstep(h,t,z,T,approx);
		moultonstep(h,t,z,T,approx);
		
		fprintf(out,"%lf %lf %lf\n",t[4],z[2][4],z[0][4]);
	
		if(fabs(z[2][4]-x0)/h < 1 && fabs(z[0][4] - y0)/h < 1 && flag == 1)
		{
			periodo = t[4];
			flag = 0;
		}
		
		r = pow(z[2][4]*z[2][4] + z[0][4]*z[0][4],1.0/2.0);
		raverage += r;
		if(r > da)
		{
			da = r;
		}
		if(r < dp)
		{

			dp = r;
		}
		
		
		
		for(i = 0;i < 5;i++)
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
	stats[0] = (da - dp)/(da + dp);
	stats[1] = periodo;
	raverage = raverage/avg;
	rerror = fabs(raverage - da);
	if(rerror > fabs(raverage - dp))
	{
		rerror = fabs(raverage - dp);
	}
	printf("r average = %lf +- %lf\n",raverage,rerror);
	
	fclose(out);
	return 0.0;	
}
void rk4(double h,double t[],double z[][5])
{
	double k[4];
	double vect[4];
	double z_temp[4][5] = {};
	int i,j,r,ind;
	for(j = 1; j < 4;j++)
	{
		ind = j-1;
		for(i = 0;i < 4;i++)
		{
			g(t[j-1],z,vect,ind);
			k[0] = vect[i];
			
			for(r = 0;r < 4;r++)
			{
				z_temp[r][ind] = z[r][ind] + 0.5*h*k[0];
			} 
			
			g(t[j-1] + 0.5*h,z_temp,vect,ind);
			k[1] = vect[i];
			for(r = 0;r < 4;r++)
			{
				z_temp[r][ind] = z[r][ind] + 0.5*h*k[1];
			} 
			
			g(t[j-1] + 0.5*h,z_temp,vect,ind);
			k[2] = vect[i];
			
			for(r = 0;r < 4;r++)
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
//signle step of the bashcroft method
	int i;
	double vect0[4] = {},vect1[4] = {},vect2[4] = {},vect3[4] = {};
	g(t[0],z,vect0,0);
	g(t[1],z,vect1,1);
	g(t[2],z,vect2,2);
	g(t[3],z,vect3,3);

	for(i = 0;i < 4;i++)
	{
		z[i][4] = z[i][3] + h/24.0*(55*vect3[i] - 59*vect2[i] + 37*vect1[i] - 9*vect0[i]);
		approx[i][0] = z[i][4];
	}
	
	h = h/2;
	for(i = 0;i < 4;i++)
	{
		approx[i][1] = z[i][3] + h/24.0*(55*vect3[i] - 59*vect2[i] + 37*vect1[i] - 9*vect0[i]);
	}
	
	
}
void moultonstep(double h,double t[],double z[][5],double T,double approx[][2])
{
//single step of the multon method
    int i;
    double x0,x1,x2,f0,f1;
    double vect1[4] = {},vect2[4] = {},vect3[4] = {},vect4[4] = {},vect41[4] = {};
    g(t[1],z,vect1,1);
    g(t[2],z,vect2,2);
    g(t[3],z,vect3,3);
    g(t[4],z,vect4,4);
     
    for(i = 0;i < 4;i++)
    {
        z[i][4] = approx[i][1];
    }
    g(t[4],z,vect41,4);
     
     
    for(i = 0;i < 4;i++)
    {
        x0 = approx[i][0];
        x1 = approx[i][1];
         
        f0 = approx[i][0] - z[i][3] - h/24.0*(9*vect4[i] + 19*vect3[i] - 5*vect2[i] + vect1[i]);
        f1 = approx[i][1] - z[i][3] - h/24.0*(9*vect41[i] + 19*vect3[i] - 5*vect2[i] + vect1[i]);
         
        x2 = x1 - (x1-x0)/(f1-f0)*f1;
        z[i][4] = x2;
    }
     
}
void g(double t,double z[][5],double vect[],int ind)
{
	double y = z[0][ind];
	double vy = z[1][ind];
	double x = z[2][ind];
	double vx = z[3][ind];

	vect[0] = vy;
	vect[1] = -y/pow(x*x + y*y,1.5);
	vect[2] = vx;
	vect[3] = -x/pow(x*x + y*y,1.5);
}
