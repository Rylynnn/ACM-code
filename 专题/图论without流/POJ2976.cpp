#include <iostream>
#include <cstring>
#include <stack>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <queue>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <deque>
#include <bitset>
#include <algorithm>
  
using namespace std;
  
const double eps=1e-10;
const double pi=3.1415926535897932384626433832795;
const double eln=2.718281828459045235360287471352;
  
#define LL long long
#define IN freopen("in.txt", "r", stdin)
#define OUT freopen("out.txt", "w", stdout)
#define scan(x) scanf("%d", &x)
#define scan2(x, y) scanf("%d%d", &x, &y)
#define scan3(x, y, z) scanf("%d%d%d", &x, &y, &z)
#define sqr(x) (x) * (x)
#define pr(x) printf("Case %d: ",x)
#define prn(x) printf("Case %d:\n",x)
#define prr(x) printf("Case #%d: ",x)
#define prrn(x) printf("Case #%d:\n",x)
#define lc o << 1
#define rc o << 1 | 1
#define lowbit(x) (x&(-x))
#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))

const int maxn=1005;
double a[maxn],b[maxn];
struct hei
{
	double num;
	int pos;
	hei(double num=0,int pos=0):num(num),pos(pos){}
	bool operator < (struct hei p)const
	{return num>p.num;}
}d[maxn];
double p,q,ans,l;
int i,n,m;

int main()
{
	while(scan2(n,m)==2 && n+m>0)
	{
		m=n-m;
		for(i=1;i<=n;i++)scanf("%lf",&a[i]);
		for(i=1;i<=n;i++)scanf("%lf",&b[i]);
		l=0;
		while(true)
		{
			ans=l;
			for(i=1;i<=n;i++)d[i]=hei(a[i]-ans*b[i],i);
			sort(d+1,d+n+1);
			double fz,fm;
			fz=fm=0.0;
			for(i=1;i<=m;i++)
			{
				fz+=a[d[i].pos];
				fm+=b[d[i].pos];
			}
			l=fz/fm;
			if(fabs(ans-l)<eps)break;
		}
		printf("%.0f\n",100.0*ans);
	}
    return 0;
}

