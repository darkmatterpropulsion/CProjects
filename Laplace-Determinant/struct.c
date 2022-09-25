#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

struct matrix
{
  int nLines;
  int nColumns;
  int ** matrice;
};


void fillMatrix(struct matrix* A);
void createMatrix(struct matrix* A);
long long int calculateDeterminant(struct matrix A);
void fillMinor(struct matrix *minor,struct matrix A,int removedLine,int removedColumn);
void printMatrix(struct matrix* A);


int main()
{
  struct matrix A;
  
  printf("Insert the number of lines for matrix1: ");
  scanf("%d",&A.nLines);
  printf("\nInsert the number of columns for matrix1:  ");
  scanf("%d",&A.nColumns);
  


  createMatrix(&A);  //  I create the 3 matrix;   

  fillMatrix(&A);   //   I fill 2 of these. matrix C is used to store sum value.
  printMatrix(&A);
  printf("\n\n\n");
 

  long long int det = calculateDeterminant(A);
  printf("\n the Determinant is : %lld \n",det);
  
  
  for(int i = 0; i < A.nLines; i++)
  	free(A.matrice[i]);
  free(A.matrice);


  return 0;

}

void createMatrix(struct matrix *A)
{
  //create the matrix
  A->matrice = (int **)calloc(A->nLines,sizeof(int*));  
  for(int i=0; i < A->nLines; i++)
    {
      A->matrice[i] = (int *)calloc(A->nColumns,sizeof(int));
    }
}


void fillMatrix(struct matrix * A) //this routine serves to fill the matrix with some int
				   //so that we can readily see the program work. 
{
  srand(time(NULL));
  int i,j;
  for(i=0;i<A->nLines;i++)
    for(j=0;j<A->nColumns;j++)
      {
	A->matrice[i][j] = rand()%99;
      }
}

long long int calculateDeterminant(struct matrix A)
{
//This is the main recursive function that gives the determinant.
  if(A.nLines != A.nColumns)
    {
      //The determinant is defined only for a square matrix.
      printf("Not square matrix.");
      return -1;
    }
    
  int j=0;
  long long int temp=0; //Will store the partial sum of the determinant for each step.
  struct matrix minor;
  for(j=0;j < A.nColumns;j++) //the laplace determinant is calculated by column, on the first line.
  {
	  if(A.nColumns==2 && A.nLines==2) //The trivial case of the recursive call is the 2x2 matrix
	   {
	     	 return (A.matrice[0][0]*A.matrice[1][1] - A.matrice[0][1]*A.matrice[1][0]);
	   }
	  
	  minor.nLines = A.nLines-1;       //setting up the dimension of the minor.
	  minor.nColumns = A.nColumns-1;
	  createMatrix(&minor); //allocating the memory for the minor.
	  fillMinor(&minor,A,0,j); //filling the minor with the remaining elements of the original matrix.
	  printf("\n\n");
	  printMatrix(&minor); // it print the minor, given that it's the main focus of this exercise.
	  
	  //the actual recursive call.
	  if((2+j) % 2 == 0)
	  {
	  	temp += A.matrice[0][j]*calculateDeterminant(minor);
	  }else
	  {
	  	temp += -1*A.matrice[0][j]*calculateDeterminant(minor);
	  }
  }
  
  for(int i = 0; i < minor.nLines; i++)
  	free(minor.matrice[i]);
  free(minor.matrice);
  
  return temp;
  
  
}

void fillMinor(struct matrix * minor,struct matrix A,int removedLine,int removedColumn)
{
//for a removed line and a removed column we calculate the minor using the plain
//mathematical definition.

  int i,j;
  for(i=0;i<A.nLines;i++)
    {
      if(i == removedLine)
	continue;
      
      for(j=0;j<A.nColumns;j++)
	{
	  if(j < removedColumn)
	    minor->matrice[i-1][j] = A.matrice[i][j];
	  else if(j == removedColumn)
	    continue;
	  else if(j > removedColumn)
	    minor->matrice[i-1][j-1]=A.matrice[i][j];
	} 
    }
}


void printMatrix(struct matrix *A)
{
//service routine to print the matrix
  int h,k;
  for(h=0;h<A->nLines;h++)
    {
      for(k=0;k<A->nColumns;k++)
	{
	  printf(" %d ",A->matrice[h][k]);
	}
      printf("\n");
    }
}
