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
#define other(x,y) ((edg[x].p==y)?(edg[x].q):(edg[x].p))

const int maxn=20005;
const int maxm=100005;

struct edge
{
	int p,q;
	edge(int p=0,int q=0):p(p),q(q){}
}edg[maxm];

vector<int> g[maxn];
int bcc[maxm],pre[maxn];
int p[maxm],s[maxm];
int dfs_clock,bcc_cnt;

int dfs(int u,int fa)
{
	int lowu=pre[u]=++dfs_clock;
	for(unsigned int i=0;i<g[u].size();i++)
	{
		int side=g[u][i];
		int v=other(side,u);
		if(!pre[v])
		{
			s[++s[0]]=side;
			int lowv=dfs(v,u);
			lowu=min(lowu,lowv);
			if(lowv>=pre[u])
			{
				bcc_cnt++;
				for(;;)
				{
					int x=s[s[0]];s[0]--;
					bcc[x]=bcc_cnt;
					p[bcc_cnt]=min(p[bcc_cnt],x);
					if(x==side)break;
				}
			}
		}else if(pre[v]<pre[u] && v!=fa)
		{
			s[++s[0]]=side;lowu=min(lowu,pre[v]);
		}
	}
	return lowu;
}

void tarjan(int n)
{
	dfs_clock=bcc_cnt=0;
	memset(pre,0,sizeof(pre));
	memset(bcc,0,sizeof(bcc));
	memset(p,0x3f3f3f3f,sizeof(p));
	s[0]=0;
	for(int i=1;i<=n;i++)if(!pre[i])dfs(i,-1);
}

int main()
{
	int n,m;
	scan2(n,m);
	for(int i=1;i<=m;i++)
	{
		int j,k;
		scan2(j,k);
		edg[i]=edge(j,k);
		g[j].push_back(i);
		g[k].push_back(i);
	}
	tarjan(n);
	printf("%d\n",bcc_cnt);
	for(int i=1;i<=m;i++)
	{
		if(i>1)printf(" ");
		printf("%d",p[bcc[i]]);
	}
    return 0;
}

