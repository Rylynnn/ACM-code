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
#define type double

const type inf=2000000000.0;
const int maxn=105;

struct edge
{
	int from,to;
	type cost;
	edge(int from=0,int to=0,double cost=0):from(from),to(to),cost(cost){}
}edg[10005];

struct point
{
	type p,q;
	point(type p=0,type q=0):p(p),q(q){}
}a[maxn];

int pre[maxn],id[maxn],vis[maxn];
type in[maxn];

type ZLEdmonds(int n,int m,int root)//自环在输入建图时直接忽略，如需加入，可另存
{
	type tot=0.0;
	//判断是否有树
	while(true)
	{
		for(int i=1;i<=n;i++)in[i]=inf;
		for(int i=1;i<=m;i++)
		{
			int u=edg[i].from;
			int v=edg[i].to;
			if(edg[i].cost<in[v] && u!=v){pre[v]=u;in[v]=edg[i].cost;}
		}
		for(int i=1;i<=n;i++)if(i!=root && in[i]==inf)return -1;
		//找环
		int cnt=1;
		memset(id,0,sizeof(id));
		memset(vis,0,sizeof(vis));
		in[root]=0;
		for(int i=1;i<=n;i++)//标记每个环
		{
			tot+=in[i];
			int v=i;
			while(vis[v]!=i && id[v]==0 && v!=root)
			{vis[v]=i;v=pre[v];}
			if(v!=root && id[v]==0)//缩点
			{
				for(int u=pre[v];u!=v;u=pre[u])id[u]=cnt;
				id[v]=cnt++;
			}
		}
		if(cnt==1)break;
		for(int i=1;i<=n;i++)if(id[i]==0)id[i]=cnt++;
		//建立新图
		for(int i=1;i<=m;i++)
		{
			int u=edg[i].from;
			int v=edg[i].to;
			edg[i].from=id[u];
			edg[i].to=id[v];
			if(id[u]!=id[v])edg[i].cost-=in[v];
		}
		n=cnt-1;
		root=id[root];
	}
	return tot;
}

type js(point x,point y)
{
	return sqrt((x.p-y.p)*(x.p-y.p)+(x.q-y.q)*(x.q-y.q));
}

int main()
{
	int i,j,k,l,m,n;
	while(scan2(n,m)!=EOF)
	{
		memset(a,0,sizeof(a));
		memset(edg,0,sizeof(edg));
		l=0;
		for(i=1;i<=n;i++)scanf("%lf %lf",&a[i].p,&a[i].q);
		for(i=1;i<=m;i++)
		{
			scan2(j,k);
			if(j==k)continue;
			edg[++l]=edge(j,k,js(a[j],a[k]));
		}
		type ans=ZLEdmonds(n,l,1);
		if(ans<0)printf("poor snoopy\n");else printf("%.2f\n",ans);
	}
    return 0;
}

