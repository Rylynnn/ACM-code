#include<bits/stdc++.h>

const int maxn=1005;
const double eps=1e-8;
#define LL long long int

using namespace std;

int phi[maxn];

void swap(double& p,double& q)
{
	double t;
	t=p;p=q;q=t;
}

struct Matrix
{
	double a[maxn][maxn];
	//1-n�б�ʾ��1-n������
	//ÿ�е�1-n��Ԫ�ر�ʾϵ������n+1��Ԫ�ر�ʾ�Ⱥ��ұߵĳ���
}q;

int ii,jj,nn;

LL det(LL a[][maxn], int n) {//������ʽֵ(������)
    int i, j, k, r;
    LL res = 1;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            while (a[j][i]) {
                LL f = a[i][i] / a[j][i];
                for (int k = i; k < n; k++) a[i][k] -= f * a[j][k];
                for (int k = i; k < n; k++) swap(a[i][k], a[j][k]);
                res = -res;
            }
        }
        if (a[i][i] == 0) return 0;
        res *= a[i][i];
    }
    return res < 0 ? -res : res;
}

double FF(double x)//����ֵĺ����������޸�
{
	return 1.0;
}

double simpson(double x,double y)
{
	double z=x+(y-x)/2.0;
	return (y-x)/6.0*(FF(x)+FF(y)+4*FF(z));
}

double asr(double x,double y,double eeps,double A)//eepsΪ����
{
	double z=x+(y-x)/2.0;
	double L=simpson(x,z);
	double R=simpson(z,y);
	if(fabs(L+R-A)<=15*eeps)return (L+R)+(L+R-A)/15.0;
	else return asr(x,z,eeps/2.0,L)+asr(z,y,eeps/2.0,R);
}

double simpson_zsx(double x,double y,double eeps)//����Ӧ����ɭ������
{
	return asr(x,y,eeps,simpson(x,y));
}

void gauss_eli(struct Matrix& p,int n)//��˹��Ԫ
{
	int i,j,k,r;
	for(i=1;i<=n;i++)
	{
		r=i;
		for(j=i+1;j<=n;j++)
			if(fabs(p.a[j][i])>fabs(p.a[r][i]))r=j;
		if(r!=i)for(j=1;j<=n+1;j++)swap(p.a[r][j],p.a[i][j]);
		for(k=1;k<=i-1;k++)
		{
			if(p.a[i][k]==0)continue;
			for(j=n+1;j>=k;j--)
				p.a[i][j]-=p.a[k][j]/p.a[k][k]*p.a[i][k];
		}
	}
	for(i=n;i>=1;i--)
	{
		for(j=i+1;j<=n;j++)
			p.a[i][n+1]-=p.a[j][n+1]*p.a[i][j];
		p.a[i][n+1]/=p.a[i][i];
	}
}

LL gcd(LL a,LL b)
{
	return b==0?a:gcd(b,a%b);
}

void tgcd(LL a,LL b,LL& d,LL& x,LL& y)//��չŷ�����
{
	if(!b){d=a;x=1;y=0;}
	else{tgcd(b,a%b,d,y,x);y-=x*(a/b);}
}

LL pow_mod(LL a,LL p,LL n)//ͬ�������
{
	if(p==0)return 1;
	LL ans=pow_mod(a,p/2,n);
	ans=(ans*ans)%n;
	if(p%2==1)ans=(ans*a)%n;
	return ans;
}

int euler_phi(int n)//��ŷ������
{
	int m=(int)sqrt(n+0.5);
	int ans=n;
	for(int i=2;i<=m;i++)
		if(n%i==0)
		{
			ans=ans/i*(i-1);
			while(n%i==0)n=n/i;
		}
	if(n>1)ans=ans/n*(n-1);
	return ans;
}

void phi_table(int n)//ŷ��������
{
	memset(phi,0,n+1);
	phi[1]=1;
	for(int i=2;i<=n;i++)
	{
		if(phi[i])continue;
		for(int j=i;j<=n;j+=i)
		{
			if(!phi[j])phi[j]=j;
			phi[j]=phi[j]/i*(i-1);
		}
	}
}

LL inv(LL a,LL n)//a����n����Ԫ
{
	LL d,x,y;
	tgcd(a,n,d,x,y);
	return d==1?(x+n)%n:-1;
}

LL china(int n,int* a,int* m)//�й�ʣ�ඨ��
{
	LL M=1,d,y,x=0;
	for(int i=0;i<n;i++)M*=m[i];
	for(int i=0;i<n;i++)
	{
		LL w=M/m[i];
		tgcd(m[i],w,d,d,y);
		x=(x+y*w*a[i])%M;
	}
	return (x+M)%M;
}

int log_mod(int a,int b,int n)//���ģ����a^x=b(mod n),nΪ����,�޽ⷵ��-1
{
	int m,v,e=1,i;
	m=(int)sqrt(n+0.5);
	v=inv(pow_mod(a,m,n),n);
	map<int,int> x;
	x[1]=0;
	for(i=1;i<m;i++)
	{
		e=(e*a)%n;
		if(!x.count(e))x[e]=i;
	}
	for(i=0;i<m;i++)
	{
		if(x.count(b))return (i*m+x[b]);
		b=(b*v)%n;
	}
	return -1;
}



int main()
{
	freopen("in.txt","r",stdin);
	scanf("%d",&nn);
	for(ii=1;ii<=nn;ii++)
		for(jj=1;jj<=nn+1;jj++)scanf("%lf",&q.a[ii][jj]);
	gauss_eli(q,nn);
	for(ii=1;ii<=nn;ii++)
	{
		printf("a%d=%0.5lf",ii,q.a[ii][nn+1]);
		printf("\n");
	}
	return 0;
}