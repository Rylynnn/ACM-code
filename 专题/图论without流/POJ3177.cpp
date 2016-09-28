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

const int maxn=5005;

struct edge
{
	int p,q;
	edge(int p=0,int q=0):p(p),q(q){}
}edg[2*maxn];
bool bridge[2*maxn];
vector<int> g[maxn];
int pre[maxn],dfs_clock,bcc_cnt;
int i,j,k,l,m,n;

int other(int t,int x)
{
	return (edg[t].p==x?edg[t].q:edg[t].p);
}

int dfs(int u,int fa)
{
	int lowu=pre[u]=++dfs_clock;
	for(unsigned int i=0;i<g[u].size();i++)
	{
		int v=other(g[u][i],u);
		if(!pre[v])
		{
			int lowv=dfs(v,u);
			if(lowv>pre[u])bridge[g[u][i]]=true;
			lowu=min(lowu,lowv);
		}else if(pre[v]<pre[u] && v!=fa)lowu=min(lowu,pre[v]);
	}
	return lowu;
}

void find_bcc(int x)
{
	pre[x]=bcc_cnt;
	for(unsigned int i=0;i<g[x].size();i++)
	{
		int side=g[x][i];
		if(bridge[side])continue;
		int v=other(side,x);
		if(!pre[v])find_bcc(v);
	}
}

void tarjan(int n)
{
	dfs_clock=bcc_cnt=0;
	memset(pre,0,sizeof(pre));
	for(int i=1;i<=n;i++)if(!pre[i])dfs(i,-1);
	memset(pre,0,sizeof(pre));
	for(int i=1;i<=n;i++)
	{
		if(pre[i])continue;
		bcc_cnt++;
		find_bcc(i);
	}
}

int main()
{
	scan2(n,m);
	for(i=1;i<=m;i++)
	{
		scan2(j,k);
		edge e=edge(j,k);
		edg[i]=e;
		g[j].push_back(i);
		g[k].push_back(i);
	}
	memset(bridge,0,sizeof(bridge));
	tarjan(n);
	int tot=0;
	int sum[bcc_cnt+1];
	bool ss[bcc_cnt+1][bcc_cnt+1];
	memset(ss,0,sizeof(ss));
	memset(sum,0,sizeof(sum));
	for(i=1;i<=m;i++)
	{
		int pp=edg[i].p;
		int qq=edg[i].q;
		pp=pre[pp],qq=pre[qq];
		if(pp==qq)continue;
		if(!ss[pp][qq]){sum[pp]++;sum[qq]++;}
		ss[pp][qq]=ss[qq][pp]=true;
	}
	for(i=1;i<=bcc_cnt;i++)if(sum[i]==1)tot++;
	printf("%d",(tot+1)/2);
    return 0;
}

