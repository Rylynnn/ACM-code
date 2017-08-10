/*
    Author:Scarlet
*/
#include<bits/stdc++.h>
#define maxn 262144
using namespace std;
typedef long long LL;
#define G c=getchar()
inline int read()
{
    int x=0,f=1;char G;
    while(c>57||c<48){if(c=='-')f=-1;G;}
    while(c>47&&c<58){x=x*10+c-48;G;}
    return x*f;
}
#define pi M_PI
typedef complex<double>C;
int N,g[maxn];
C a[maxn],b[maxn];
void DFT(C *a,int f)
{
    for(int i=0;i<N;i++)if(g[i]>i)swap(a[i],a[g[i]]);
    for(int i=1;i<N;i<<=1)
    {
        C e(cos(pi/i),f*sin(pi/i));
        for(int j=0;j<N;j+=i<<1)
        {
            C w(1,0);
            for(int k=0;k<i;k++,w*=e)
            {
                C x=a[j+k],y=w*a[j+k+i];
                a[j+k]=x+y;a[j+k+i]=x-y;
            }
        }
    }
    if(f-1)for(int i=0;i<N;i++)a[i]/=N;
}
void mul(C *a,C *b,int n)
{
    int t=-1;
    for(N=1;N<=n;N<<=1,t++);
    for(int i=1;i<N;i++)g[i]=(g[i>>1]>>1)|((i&1)<<t);
    DFT(a,1);DFT(b,1);
    for(int i=0;i<N;i++)
        a[i]=a[i]*b[i];
    DFT(a,-1);
}
int main()
{
    int n=read(),m=read();
    for(int i=0;i<=n;i++)a[i]=read();
    for(int i=0;i<=m;i++)b[i]=read();
    mul(a,b,n+m+1);
    for(int i=0;i<=n+m;i++)
        printf("%d ",(int)(a[i].real()+0.5));
}
