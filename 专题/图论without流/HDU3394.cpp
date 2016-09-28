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

struct edge
{
	int p,q;
	edge(int p=0,int q=0):p(p),q(q){}
};

vector<int> g[maxn];
vector<int> bcc[maxn];
int bccno[maxn];
unsigned int side[maxn];
int pre[maxn],dfs_clock,bcc_cnt;
int i,j,k,l,m,n;
int ans;
stack<edge> s;

int dfs(int u,int fa)
{
	int lowu=pre[u]=++dfs_clock;
	int child=0;
	edge e;
	for(unsigned int i=0;i<g[u].size();i++)
	{
		int v=g[u][i];
		e=edge(u,v);
		if(!pre[v])
		{
			s.push(e);
			child++;
			int lowv=dfs(v,u);
			lowu=min(lowv,lowu);
			if(lowv>=pre[u])
			{
				if(lowv>pre[u])ans++;
				bcc_cnt++;bcc[bcc_cnt].clear();
				for(;;)
				{
					edge x=s.top();s.pop();
					side[bcc_cnt]++;
					if(bccno[x.p]!=bcc_cnt){bcc[bcc_cnt].push_back(x.p);bccno[x.p]=bcc_cnt;}
					if(bccno[x.q]!=bcc_cnt){bcc[bcc_cnt].push_back(x.q);bccno[x.q]=bcc_cnt;}
					if(x.p==u && x.q==v)break;
				}
			}
		}else if(pre[v]<pre[u] && v!=fa)
		{
			s.push(e);lowu=min(lowu,pre[v]);
		}
	}
	return lowu;
}

void tarjan(int n)
{
	memset(bccno,0,sizeof(bccno));
	memset(pre,0,sizeof(pre));
	memset(side,0,sizeof(side));
	dfs_clock=bcc_cnt=ans=0;
	for(i=0;i<n;i++)if(!pre[i])dfs(i,-1);	
}

int main()
{
	while(scan2(n,m)==2 && n+m>0)
	{
		for(i=0;i<=n;i++)g[i].clear();
		for(i=1;i<=m;i++)
		{
			scan2(j,k);
			g[j].push_back(k);
			g[k].push_back(j);
		}
		tarjan(n);
		int tot=0;
		for(i=1;i<=bcc_cnt;i++)
			if(side[i]>bcc[i].size())tot+=side[i];
		printf("%d %d\n",ans,tot);
	}
    return 0;
}

