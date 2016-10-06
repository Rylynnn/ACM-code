#include<bits/stdc++.h>

using namespace std;

const double eps=1e-10;
const double pi=3.1415926535897932384626433832795;
const double eln=2.718281828459045235360287471352;

const int maxn=105000;

complex<double> epsilon[maxn];
complex<double> arti_epsilon[maxn];
complex<double> a[maxn],b[maxn],c[maxn],temp[maxn];

int n1,n2,m;

void init_epsilon(int n)
{
	for(int i=0;i!=n;i++)
	{
		epsilon[i]=complex<double>(cos(2.0*pi*i/n),sin(2.0*pi*i/n));
		arti_epsilon[i]=conj(epsilon[i]);
	}
}

int calc(int t)
{
	int j=0;
	while((1<<j)<=t)j++;
	return 1<<j;
}


void DFT(int n,complex<double>*  buffer,int offset,int step,complex<double>* epsilon)
{
	if(n==1)return;
	int m=n>>1;
	DFT(m,buffer,offset,step<<1,epsilon);
	DFT(m,buffer,offset+step,step<<1,epsilon);
	for(int k=0;k!=m;k++)
	{
		int pos=2*step*k;
		temp[k]=buffer[pos+offset]+epsilon[k*step]*buffer[pos+offset+step];
		temp[k+m]=buffer[pos+offset]-epsilon[k*step]*buffer[pos+offset+step];
	}
	for(int i=0;i!=n;i++)buffer[i*step+offset]=temp[i];
}

//IDFT将DFT的epsilon改为arti_epsilon

void FFT(int m,complex<double>* a,complex<double>* b,complex<double>* c)
{
	init_epsilon(m);
	DFT(m,a,0,1,epsilon);
	DFT(m,b,0,1,epsilon);
	for(int i=0;i<=m;i++)c[i]=a[i]*b[i];
	IDFT(m,c,0,1,epsilon);
	double mm=m;
	for(int i=0;i<=m;i++)c[i]/=mm;
}

int init()//n1,n2表示多项式次数
{
	double x,y;
	scanf("%d%d",&n1,&n2);
	memset(a,0,sizeof(a));
	memset(b,0,sizeof(b));
	for(int i=0;i<=n1;i++)
	{
		scanf("%lf %lf",&x,&y);
		a[i].real(x);
		a[i].imag(y);
	}
	for(int i=0;i<=n2;i++)
	{
		scanf("%lf %lf",&x,&y);
		b[i].real(x);
		b[i].imag(y);
	}
	m=calc(n1+n2);
	return m;
}

void print()
{
	for(int i=0;i<m;i++)printf("%lf %lf\n",real(c[i]),imag(c[i]));
}

int main()
{
	m=init();
	FFT(m,a,b,c);
	print();
}
