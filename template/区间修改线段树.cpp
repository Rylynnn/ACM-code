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

const int maxn=1000005;

struct tree
{
	int mi,ls,rs,ll,rr,add;
}a[3*maxn];

int c[maxn];
int i,j,k,l,m,n,t,T;

void update(int x)
{
	if(a[x].ls+a[x].rs==0)return;
	a[x].mi=min(a[a[x].ls].mi,a[a[x].rs].mi);
}

void downdate(int x)
{
	if(a[x].ls+a[x].rs==0)return;
	if(a[x].add==0)return;
	a[a[x].ls].add+=a[x].add;
	a[a[x].ls].mi+=a[x].add;
	a[a[x].rs].add+=a[x].add;
	a[a[x].rs].mi+=a[x].add;
	a[x].add=0;
}

void mt(int l,int r)
{
	if(l==r)
	{
		a[k].ll=a[k].rr=l;
		a[k].mi=c[l];
		a[k].ls=a[k].rs=0;
		return;
	}
	int t=k;
	int mid=(l+r)>>1;
	a[t].ll=l;a[t].rr=r;
	k++;a[t].ls=k;mt(l,mid);
	k++;a[t].rs=k;mt(mid+1,r);
	update(t);
}

void add(int l,int r,int nu,int x)
{
	if(a[x].ll==l && a[x].rr==r)
	{
		a[x].add+=nu;
		a[x].mi+=nu;
		return;
	}
	downdate(x);
	int mid=(a[x].ll+a[x].rr)>>1;
	if(mid<l)add(l,r,nu,a[x].rs);
	else if(mid>=r)add(l,r,nu,a[x].ls);
	else
	{
		add(l,mid,nu,a[x].ls);
		add(mid+1,r,nu,a[x].rs);
	}
	update(x);
}

int find(int l,int r,int x)
{
	if(a[x].ll==l && a[x].rr==r)return a[x].mi;
	downdate(x);
	update(x);
	int mid=(a[x].ll+a[x].rr)>>1;
	if(mid<l)return find(l,r,a[x].rs);
	if(mid>=r)return find(l,r,a[x].ls);
	return min(find(l,mid,a[x].ls),find(mid+1,r,a[x].rs));
}

int main()
{
	scanf("%d%d",&n,&m);
	for(i=1;i<=n;i++)scanf("%d",&c[i]);
	k=1;mt(1,n);
	for(i=1;i<=m;i++)
	{
		scan3(l,j,k);
		if(find(j,k,1)<l)
		{
			printf("-1\n%d",i);
			break;
		}
		add(j,k,-l,1);
	}
	if(i>m)printf("0");
	return 0;
}


