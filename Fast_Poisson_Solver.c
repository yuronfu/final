#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int Print_Matrix(double **A, int N)
{
	int i, j;
	for(i=0;i<N;++i)
	{
		for(j=0;j<N;++j)
		{
			printf("%.5f ",A[i][j]);
		}
		printf("\n");
	}
}
int Exact_Solution(double **U, int N)
{
	// put the exact solution 
	int i,j,k;
	double x, y, h;
	h = 1.0/N;
	for(i=0;i<N-1;++i)
	{
		for(j=0;j<N-1;++j)
		{
			//k = j + i*(N-1);
			x = (i+1)*h;
			y = (j+1)*h;
			U[i][j] = sin(M_PI*x)*sin(2*M_PI*y);
		}
	}
	return 1;
}
int Exact_Source(double **F, int N)
{
	int i,j,k;
	double x, y, h;
	h = 1.0/N;
	for(i=0;i<N-1;++i)
	{
		for(j=0;j<N-1;++j)
		{
			//k = j + i*(N-1);
			x = (i+1)*h;
			y = (j+1)*h;
			F[i][j] = -(1.0+4.0)*h*h*M_PI*M_PI*sin(M_PI*x)*sin(2*M_PI*y);
		}
	}
	return 1;	
}
double Error(double *x, double *u, int N)
{
	// return max_i |x[i] - u[i]|
	int i;
	double v = 0.0, e;
	
	for(i=0;i<N;++i)
	{
		e = fabs(x[i]-u[i]);
		if(e>v) v = e;
	}
	return v;
}
void bit_reverse(double *x_r,double *x_i,double *y_r,double *y_i,int N)
{
     int i,j,M,t;
     
     for(i = 0 ; i < N ; i++)
     {
           t = i;
           M = N/2;
           j = 0;
           while(M >=1)
           {
                   if(t/M > 0)
                   {
                          j = j + (t/M)*N/(2*M);
                          t = t - M;
                   }
                   M = M/2;
           }
           //y[i] = x[j]
           y_r[i] = x_r[j];
           y_i[i] = x_i[j];
     }
     
     return;
}
void FFT(double *x_r,double *x_i,double *y_r,double *y_i,int N)
{
     int i,j,t,p = 2;
     double theta,w_r,w_i,temp_r,temp_i;
     
     bit_reverse(x_r,x_i,y_r,y_i,N);
     
     for(i = 1 ; i < N ; i*=p)
     {
           t = 0;
           
           theta = -M_PI/i;
           while(t < N)
           {
                   for(j = 0 ; j < i ; j++)
                   {
                         w_r = cos(j*theta);
                         w_i = sin(j*theta);
                         temp_r = w_r*y_r[j+i+t] - w_i*y_i[j+i+t];
                         temp_i = w_i*y_r[j+i+t] + w_r*y_i[j+i+t];
                         y_r[j+i+t] = y_r[j+t] - temp_r;
                         y_i[j+i+t] = y_i[j+t] - temp_i;
                         y_r[j+t]+=temp_r;
                         y_i[j+t]+=temp_i;
                   }
                   t = t + p*i;
           }
     }
     
     
     return;
}
void DST(double *x,double *y,int N)
{
     double *x_r,*x_i,*y_r,*y_i;
     int k;
     
     x_r = (double *) malloc (N*sizeof(double));
     x_i = (double *) malloc (N*sizeof(double));
     
     y_r = (double *) malloc (N*sizeof(double));
     y_i = (double *) malloc (N*sizeof(double));
     
     for(k = 0 ; k < N ; k++)
     {
           *(x_r+k) = 0;
           *(x_i+k) = 0;
     }
     for(k = 0 ; k < N/2-1 ; k++)
     {
           x_r[k+1] = x[k];
           x_r[k+N/2+1] = -x[N/2-2-k];
     }
     //for(k = 0 ; k < N ; k++) printf("x_%d : %f + %f i\n",k,x_r[k],x_i[k]);
     FFT(x_r,x_i,y_r,y_i,N);
     //for(k = 0 ; k < N ; k++) printf("y_%d : %f + %f i\n",k,y_r[k],y_i[k]);
     
     for(k = 0 ; k < N/2-1 ; k++) y[k] = -0.5*y_i[k+1];
     
     free(x_r);
     free(x_i);
     free(y_r);
     free(y_i);
     return;
}
void iDST(double *x,double *y,int N)
{
     double *x_r,*x_i,*y_r,*y_i;
     int k;
     
     x_r = (double *) malloc (N*sizeof(double));
     x_i = (double *) malloc (N*sizeof(double));
     
     y_r = (double *) malloc (N*sizeof(double));
     y_i = (double *) malloc (N*sizeof(double));
     
     for(k = 0 ; k < N ; k++)
     {
           *(x_r+k) = 0;
           *(x_i+k) = 0;
     }
     for(k = 0 ; k < N/2-1 ; k++)
     {
           x_r[k+1] = x[k];
           x_r[k+N/2+1] = -x[N/2-2-k];
     }
     
     FFT(x_r,x_i,y_r,y_i,N);
     //for(k = 0 ; k < N ; k++) printf("y_%d : %f + %f i\n",k,y_r[k],y_i[k]);
     
     for(k = 0 ; k < N/2-1 ; k++) y[k] = -2*y_i[k+1]/N;
     
     free(x_r);
     free(x_i);
     free(y_r);
     free(y_i);
     return;
}
int Transpose(double **A, int N)
{
	int i, j;
	double v;
	for(i=0;i<N;++i)
	{
		for(j=i+1;j<N;++j)
		{
			v = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = v;
		}
	}
	return 0;
}
int DST_2D(double **X, int N)
{
	int k;
	
	for(k = 0 ; k < N ; k++) DST(X[k],X[k],2*(N+1));
	Transpose(X,N);
	for(k = 0 ; k < N ; k++) DST(X[k],X[k],2*(N+1));
    Transpose(X,N);
    
    return 0;
}
int iDST_2D(double **X, int N)
{
    int k;
	
	for(k = 0 ; k < N ; k++) iDST(X[k],X[k],2*(N+1));
	Transpose(X,N);
	for(k = 0 ; k < N ; k++) iDST(X[k],X[k],2*(N+1));
    Transpose(X,N);
    
    return 0;
}
int Fast_Poisson_Solver(double **F,int N) // N= 2^p - 1
{
    int i,j;
    double h = 1.0/(N+1);
    
    DST_2D(F,N);
    for(i = 0 ; i < N ; i++) for(j = 0 ; j < N ; j++) F[i][j] = F[i][j]/(2*(cos(M_PI*(i+1)*h)+cos(M_PI*(j+1)*h)-2));
    iDST_2D(F,N);
    
    return 0;
}
void Exact_Discretization(double **A, int N);
int GaussElimination(double **M, double *x, double *b, int N);
int main()
{
	int i, j, k, N, M; 
	double **A, *x, *u, **U, *b, **F, *r;
    
    printf("N = ");
    scanf("%d",&N);
    
    /*x = (double *) malloc (N*sizeof(double));
    y = (double *) malloc (N*sizeof(double));
    
    for(k = 0 ; k < N ; k++)
    {
          x[k] = k;
          printf("x_%d = %f\n",k,x[k]);
    }
    DST(x,y,2*(N+1));
    for(k = 0 ; k < N ; k++) printf("y_%d = %f\n",k,y[k]);
    iDST(y,y,2*(N+1));
    for(k = 0 ; k < N ; k++) printf("y_%d = %f\n",k,y[k]);
    free(x);
    free(y);*/
    
    M = (N-1)*(N-1);
    
    b = (double *) malloc(M*sizeof(double));
	// F[i][j] = F[j+i*(N-1)];
	F = (double **) malloc((N-1)*sizeof(double*));
	F[0] = b;
	for(i=1;i<N-1;++i) F[i] = F[i-1] + (N-1);
	// U[i][j] = u[j+i*(N-1)] 
	u = (double *) malloc(M*sizeof(double));
	U = (double **) malloc((N-1)*sizeof(double*));
	U[0] = u;
	for(i=1;i<N-1;++i) U[i] = U[i-1] + (N-1);
    
    /*for(i = 0 ; i < N-1 ; i++) for(j = 0 ; j < N-1 ; j++) F[i][j] = i + j;
    Print_Matrix(F, N-1);
    DST_2D(F,N-1);
    iDST_2D(F,N-1);
    Print_Matrix(F, N-1);
    system("pause");*/
    
    Exact_Source(F, N);
    
  	/*	A = (double **) malloc(M*sizeof(double*));
		A[0] = (double *) malloc(M*M*sizeof(double));
		for(i=1;i<M;++i) A[i] = A[i-1] + M;
		for(i=0;i<M*M;++i) A[0][i] = 0.0;
		x = (double *) malloc(M*sizeof(double));
		r = (double *) malloc(M*sizeof(double));
		Exact_Discretization(A,N);
		GaussElimination(A,x,b,M);
		for(i = 0 ; i < M ; i++){ printf("%.5f ",x[i]); if((i+1)%3==0) printf("\n");} */

    Exact_Solution(U, N);
    
    Fast_Poisson_Solver(F,N-1);
    printf("%e\n", Error(b, u, M));
    
    free(b);
    free(u);
    free(U);
    free(F);
    
    system("pause");
    return 0;
}
/*void Exact_Discretization(double **A, int N)
{
	int i,j,k;
	for(i=0;i<N-1;++i)
	{
		for(j=0;j<N-1;++j)
		{
			k = i + j*(N-1);
			A[k][k] = -4;
			if(j>0) A[k][k-(N-1)] = 1;
			if(i>0) A[k][k-1] = 1;
			if(i<N-2) A[k][k+1] = 1;
			if(j<N-2) A[k][k+(N-1)] = 1;
		}
	}
	return;
}
int GaussElimination(double **M, double *x, double *b, int N)
{
	// Gauss Elimination
	int i, j, k;
	double **A, v;
	
	// Copy the matrix M to A
	A = (double **) malloc(N*sizeof(double*));
	A[0] = (double *) malloc(N*N*sizeof(double));
	for(i=1;i<N;++i) A[i] = A[i-1] + N;
	for(i=0;i<N;++i)
	for(j=0;j<N;++j)
		A[i][j] = M[i][j];
	
	// put b in x 
	for(k=0;k<N;++k) x[k] = b[k];
	
	// Gauss Elimination
	for(k=0;k<N;++k)
	{
		// without pivoting
		for(i=k+1;i<N;++i)
		{
			v = A[i][k]/A[k][k];
			x[i] -= v*x[k];
			for(j=k;j<N;++j)
			{
				A[i][j] -= v*A[k][j];
			}
		}
	}
	// Backward(for upper triangular matrix) Solver
	for(i=N-1;i>=0;i--)
	{
		for(j=i+1;j<N;++j)
		{
			x[i] -= A[i][j]*x[j];
		}
		x[i] /= A[i][i];
	}
	// Free the memory of A
	free(A[0]);
	free(A);
	return;
}*/
