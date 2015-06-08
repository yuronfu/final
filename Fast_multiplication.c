#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void bit_reverse(int *x,int *y,int N)
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
void FFTn(int *x,int *y,int N)
{
     int i,j,k,t;
     long u,w,temp;
     
     bit_reverse(x,y,N);
     
     for(i = 1 ; i < N ; i*=2)
     {
           u = 1;
           if(i != 1) for(k = 0 ; k < N/(2*i) ; k++) u = (u*465)%769;
           t = 0;
           while(t < N)
           {
                   w = 1;
                   for(j = 0 ; j < i ; j++)
                   {
                         //printf("(%d,%d)\n",j+t,j+i+t);
                         //y'[j+t] = y[j+t] + w*y[j+i+t] & y'[j+i+t] = y[j+t] - w*y[j+i+t]
                         //w = w*u
                         temp = y[j+i+t]*w;
                         y[j+i+t] = (y[j+t] + 768*temp)%769;
                         y[j+t] = (y[j+t] + temp)%769;
                         w = (w*u)%769;
                   }
                   t= t + 2*i;
           }
     }
     
     return;
}
void iFFTn(int *x,int *y,int N)
{
     int i,j,k,t;
     long u,w,temp;
     
     bit_reverse(x,y,N);
     
     for(i = 1 ; i < N ; i*=2)
     {
           u = 1;
           if(i != 1) for(k = 0 ; k < N/(2*i) ; k++) u = (u*43)%769;
           t = 0;
           while(t < N)
           {
                   w = 1;
                   for(j = 0 ; j < i ; j++)
                   {
                         //printf("(%d,%d)\n",j+t,j+i+t);
                         //y'[j+t] = y[j+t] + w*y[j+i+t] & y'[j+i+t] = y[j+t] - w*y[j+i+t]
                         //w = w*u
                         temp = y[j+i+t]*w;
                         y[j+i+t] = (y[j+t] + 768*temp)%769;
                         y[j+t] = (y[j+t] + temp)%769;
                         w = (w*u)%769;
                   }
                   t= t + 2*i;
           }
     }
     
     return;
}
void split_zero_pad(long long N,int *x)
{
     int i,base = 10;
     
     for(i = 0; i < 32 ; i++)
     {
           if(N == 0) x[i] = 0;
           else
           {
               x[i] = N%base;
               N = (N - x[i])/base;
           }
           //printf("%d",x[i]);
     }
     /*
     long k=100000000,i;
     
     while(N/k == 0) k/=10;
     for(i = 0; i < 32 ; i++)
     {
           if(k == 0) x[i] = 0;
           else
           {
               x[i] = N/k;
               N = N - (N/k)*k;
               k/=10;
           }
     }
     system("pause");
     */
     return;
}
int main()
{
    int k;
    long long A,B;
    long long C;
    scanf("%lld",&A);
    scanf("%lld",&B);       /* C = A*B; printf("%lld\n",C); at most 9 digits */
    
    int *x,*y,*X,*Y;
    x = (int *) malloc (32*sizeof(int));
    y = (int *) malloc (32*sizeof(int));
    X = (int *) malloc (32*sizeof(int));
    Y = (int *) malloc (32*sizeof(int));
    split_zero_pad(A,x);
    split_zero_pad(B,y);
    
    FFTn(x,X,32);
    FFTn(y,Y,32);
    for(k = 0 ; k < 32 ; k++) X[k] = X[k]*Y[k];
    
    iFFTn(X,Y,32);
    for(k = 0 ; k < 32 ; k++)  Y[k] = (Y[k]*745)%769; //printf("%d\n",(Y[k]*745)%769);
    
    long long result = 0,base = 1;
    for(k = 0 ; k < 32 ; k++)
    {
          if(Y[k] == 0) break;
          result = result + Y[k]*base;
          base*=10;
          printf("result = %lld\n",result);
    }
    
    free(x);
    free(y);
    free(X);
    free(Y);
    system("pause");
    return 0;
}
