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
  
const int maxn=10005;
int stac[maxn],sta;
struct edge
{
	int p,q;
	edge(int p=0,int q=0):p(p),q(q){}
}edg[2*maxn];
bool used[2*maxn];
vector<int> g[maxn];
int du[maxn];
int i,j,k,l,m,n;

int other(int num,int x)
{
	return x==edg[num].p?edg[num].q:edg[num].p;
}

void dfs(int x)
{
	stac[++sta]=x;
	for(unsigned int i=0;i<g[x].size();i++)
	{
		if(used[g[x][i]])continue;
		used[g[x][i]]=true;
		dfs(other(g[x][i],x));
		break;
	}
}

void Fleury(int x)
{
	sta=1;stac[sta]=x;
	while(sta>=1)
	{
		x=stac[sta];
		bool f=false;
		for(unsigned int i=0;i<g[x].size();i++)
		{
			if(!used[g[x][i]]){f=true;break;}
		}
		if(!f)printf("%d ",stac[sta--]);
		else 
		{
			sta--;
			dfs(stac[sta+1]);
		}
	}
}

int main()
{
	//未判断图是否连通，为严谨可加入bfs检查
	scan2(n,m);
	memset(edg,0,sizeof(edg));
	for(i=1;i<=n;i++)g[i].clear();
	memset(used,0,sizeof(used));
	memset(du,0,sizeof(du));
	for(i=1;i<=m;i++)
	{
		scan2(j,k);
		edg[i]=edge(j,k);
		du[j]++;du[k]++;
		g[j].push_back(i);
		g[k].push_back(i);
	}
	int tot=0,st=0;
	for(i=1;i<=n;i++)
	{
		if(du[i]==0){tot=3;break;}
		if(du[i]%2==1){tot++;st=i;}
	}
	if(st==0)st=1;
	if(tot>2)printf("HeHeDa!");else Fleury(st);
    return 0;
}

