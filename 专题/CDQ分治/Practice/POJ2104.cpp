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

const int maxn=100005;

struct query
{
	int l,r,num;
}qu[5005];

struct hehe
{
	int pos,num;
	bool operator < (struct hehe t)const{return num<t.num;}
}a[maxn];
int i,j,k,l,m,n;
int id[maxn],c[maxn+10],tmp1[5005],tmp2[5005],ans[5005];
int mi,ma,pos=0;

void add(int x,int y)
{
	while(x<n+10)
	{
		c[x]+=y;
		x+=lowbit(x);
	}
}

int calsum(int x)
{
	int tot=0;
	while(x>0)
	{
		tot+=c[x];
		x-=lowbit(x);
	}
	return tot;
}

void solve(int h,int t,int l,int r)
{
	if(l==r || h>t)return;
	int mid=(l+r)>>1;
	while(a[pos+1].num<=mid && pos<n)
	{
		add(a[pos+1].pos,1);
		pos++;
	}
	while(a[pos].num>mid)
	{
		add(a[pos].pos,-1);
		pos--;
	}
	int l1=0,l2=0;
	for(int i=h;i<=t;i++)
	{
		if(calsum(qu[id[i]].r)-calsum(qu[id[i]].l-1)>=qu[id[i]].num)
		{
			ans[id[i]]=mid;
			tmp1[l1++]=id[i];
		}else tmp2[l2++]=id[i];
	}
	memcpy(id+h,tmp1,sizeof(int)*l1);
	memcpy(id+h+l1,tmp2,sizeof(int)*l2);
	solve(h,h+l1-1,l,mid);
	solve(h+l1,t,mid+1,r);
}

int main()
{
	scanf("%d%d",&n,&m);
	ma=0;mi=10000005;
	for(i=1;i<=n;i++)
	{
		scanf("%d",&a[i].num);
		a[i].pos=i;
		mi=min(mi,a[i].num);
		ma=max(ma,a[i].num);
	}
	sort(a+1,a+n+1);
	for(i=1;i<=m;i++)
	{
		scan3(qu[i].l,qu[i].r,qu[i].num);
		id[i]=i;
	}
	solve(1,m,mi,ma+1);
	for(i=1;i<=m;i++)printf("%d\n",ans[i]);
    return 0;
}

