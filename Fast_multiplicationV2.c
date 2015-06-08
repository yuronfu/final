#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define  MAXCHAR 100

int ReverseOrder(char *a, int N)
{
	int i, j;
	char t;
	for(i=0;i<N/2;++i)
	{
		j = N - 1 - i;
		t = a[i];
		a[i] = a[j];
		a[j] = t;
	}
	return 0;
}
int Char2Int(char *a, int N)
{
	int i;
	for(i=0;i<N;++i)
	{
		a[i] -= 48;
	}
	return 0;
}
int Int2Char(char *a, int N)
{
	int i;
	for(i=0;i<N;++i)
	{
		a[i] += 48;
	}
	a[N] = 0;
	return 0;
}
int Multiply(char *c, char *a, int Na, char *b, int Nb)
{
	int i, j;
	for(i=0;i<Na+Nb;++i) c[i] = 0;

	for(i=0;i<Na;i++)
	{
		for(j=0;j<Nb;j++)
		{
			c[i+j] += a[i]*b[j];
			if(c[i+j] >= 10)
			{
				c[i+j+1] += (c[i+j] / 10);
				c[i+j] = c[i+j] % 10;
			}
		}
	}
	if(c[Na+Nb-1]==0) return Na+Nb-1;
	else return Na+Nb;
}
int split_zero_pad(char *a, int N,long *x)
{
	int i;
	for(i=0;i<256;++i)
	{
        if(i>=N) x[i] = 0;
		else x[i] = a[i] - 48;
		//printf("%ld",x[i]);
	}
	return 0;
}
void bit_reverse(long *x,long *y,int N)
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
           y[i] = x[j];
     }
     
     return;
}
void FFTn(long *x,long *y,int N)
{
     int i,j,k,t;
     long u,w,temp;
     
     bit_reverse(x,y,N);
     
     for(i = 1 ; i < N ; i*=2)
     {
           u = 1;
           if(i != 1) for(k = 0 ; k < N/(2*i) ; k++) u = (u*7643)%9473;
           t = 0;
           while(t < N)
           {
                   w = 1;
                   for(j = 0 ; j < i ; j++)
                   {
                         //printf("(%d,%d)\n",j+t,j+i+t);
                         //y'[j+t] = y[j+t] + w*y[j+i+t] & y'[j+i+t] = y[j+t] - w*y[j+i+t]
                         //w = w*u
                         temp = (y[j+i+t]*w)%9473;
                         y[j+i+t] = (y[j+t] + 9472*temp)%9473;
                         y[j+t] = (y[j+t] + temp)%9473;
                         w = (w*u)%9473;
                   }
                   t= t + 2*i;
           }
     }
     
     return;
}
void iFFTn(long *x,long *y,int N)
{
     int i,j,k,t;
     long u,w,temp;
     
     bit_reverse(x,y,N);
     
     for(i = 1 ; i < N ; i*=2)
     {
           u = 1;
           if(i != 1) for(k = 0 ; k < N/(2*i) ; k++) u = (u*88)%9473;
           t = 0;
           while(t < N)
           {
                   w = 1;
                   for(j = 0 ; j < i ; j++)
                   {
                         //printf("(%d,%d)\n",j+t,j+i+t);
                         //y'[j+t] = y[j+t] + w*y[j+i+t] & y'[j+i+t] = y[j+t] - w*y[j+i+t]
                         //w = w*u
                         temp = (y[j+i+t]*w)%9473;
                         y[j+i+t] = (y[j+t] + 9472*temp)%9473;
                         y[j+t] = (y[j+t] + temp)%9473;
                         w = (w*u)%9473;
                   }
                   t= t + 2*i;
           }
     }
     
     return;
}
int Summation(char *c,long *x,int Na,int Nb)
{
    int i,k,n;
    
    for(i=0;i<Na+Nb;++i) c[i] = 0;
    for(k = 0; k < 256 ; k++)
    {
          if(x[k] == 0) break;
          for(n = 0 ; n < 4 ; n++)
          {
               c[k+n] += x[k]%10;
               if(c[k+n] >= 10)
               {
                         c[k+n+1] += c[k+n]/10;
                         c[k+n] = c[k+n]%10;
               }
               x[k] = x[k]/10;
          }
    }
    
   	if(c[Na+Nb-1]==0) return Na+Nb-1;
	else return Na+Nb;
}
int main()
{
	int i, j, k, Na, Nb, Nc;
	char a[MAXCHAR], b[MAXCHAR], t, c[2*MAXCHAR];
	clock_t t1, t2;
	long *x,*y,*X,*Y;
    x = (long *) malloc (256*sizeof(long));
    y = (long *) malloc (256*sizeof(long));
    X = (long *) malloc (256*sizeof(long));
    Y = (long *) malloc (256*sizeof(long));
	
	printf("Input number a and b (MAXIMUM DIGITS:%d)\n", MAXCHAR); 
	scanf("%s %s", a, b);
	printf("%s * %s = \n", a, b);
	// 算出字串a, b的長度 
	Na = strlen(a);
	Nb = strlen(b);
	// 將a, b的位元順序互換. Ex. a = 1234 -> 4321 
	ReverseOrder(a,Na);
	ReverseOrder(b,Nb);
	// 將a,b從字元值轉換成0-9的數字，可google: ascii code, 48='0', 49='1', ...  
	Char2Int(a, Na);
	Char2Int(b, Nb);
	// 做a,b的乘法(標準國小做法) 
	t1 = clock();
	Nc = Multiply(c, a, Na, b, Nb);
	t2 = clock();
	printf("%f\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	// 把c的位元順序逆序 
	ReverseOrder(c, Nc);
	// 把c的數字變回字元 
	Int2Char(c, Nc);
	printf("%s\n", c);
	
	
	Int2Char(a, Na);
	Int2Char(b, Nb);
	split_zero_pad(a, Na,x);
	split_zero_pad(b, Nb,y);
	
	FFTn(x,X,256);
    FFTn(y,Y,256);
    for(k = 0 ; k < 256 ; k++) X[k] = X[k]*Y[k];
    iFFTn(X,Y,256);
    for(k = 0 ; k < 256 ; k++) Y[k] = (Y[k]*9436)%9473; //printf("%ld\n",(Y[k]*9436)%9473);
	Nc = Summation(c,Y,Na,Nb);
	
	ReverseOrder(c, Nc); 
	Int2Char(c, Nc);
	printf("%s\n", c);
	
    free(x);
    free(y);
    free(X);
    free(Y);
	system("pause");
	return 0;
}
