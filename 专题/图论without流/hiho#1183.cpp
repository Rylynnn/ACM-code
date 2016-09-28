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

const int maxn=20005;
const int maxm=100005;

struct edge
{
	int p,q;
	edge(int p=0,int q=0):p(p),q(q){}
	bool operator < (struct edge e)const//强行相反重载
	{
		if(p>e.p)return true;
		if(p<e.p)return false;
		return q>e.q;
	}
};

vector<int> g[maxn];
int i,j,k,l,m,n;
int pre[maxn],dfs_clock;
bool iscut[maxn];

priority_queue<edge> p;

int dfs(int u,int fat)
{
	int lowu,lowv;
	lowu=pre[u]=++dfs_clock;
	int child=0;
	for(unsigned int i=0;i<g[u].size();i++)
	{
		int v=g[u][i];
		if(!pre[v])
		{
			child++;
			lowv=dfs(v,u);
			lowu=min(lowu,lowv);
			if(lowv>pre[u])p.push(edge(min(u,v),max(u,v)));
			if(lowv>=pre[u])iscut[u]=true;
		}else if(v!=fat) lowu=min(lowu,pre[v]);
	}
	if(fat==-1 && child<=1)iscut[u]=false;
	return lowu;
}

void tarjan(int n)
{
	dfs_clock=0;
	memset(pre,0,sizeof(pre));
	memset(iscut,0,sizeof(iscut));
	for(int i=1;i<=n;i++)if(!pre[i])dfs(i,-1);
}

int main()
{
	scan2(n,m);
	for(i=1;i<=m;i++)
	{
		scan2(j,k);
		g[j].push_back(k);
		g[k].push_back(j);
	}
	while(!p.empty())p.pop();
	tarjan(n);
	bool f=false;
	for(i=1;i<=n;i++)if(iscut[i]){if(f)printf(" ");printf("%d",i);f=true;}
	if(f)printf("\n");else printf("Null\n");
	while(!p.empty())
	{
		printf("%d %d\n",p.top().p,p.top().q);
		p.pop();
	}
    return 0;
}

