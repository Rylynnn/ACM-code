#include<cstdio>
#include<cstring>
#include<stack>
#include<vector>

using namespace std;

const int maxn=20005;
int pre[maxn],low[maxn],a[maxn],sccno[maxn],tot[maxn];
int edge[100005][2];
int dfs_clock,scc_cnt,maxx;
vector<int> g[maxn];
stack<int> s;

int n,m,p0,p1,i,j,k,l;

int min(int x,int y)
{
	return x<y?x:y;
}

int max(int x,int y)
{
	return x>y?x:y;
}

void dfs(int u)
{
	pre[u]=low[u]=++dfs_clock;
	s.push(u);
	for(unsigned int i=0;i<g[u].size();i++)
	{
		int v=g[u][i];
		if(!pre[v])
		{
			dfs(v);
			low[u]=min(low[u],low[v]);
		}else if(!sccno[v])low[u]=min(low[u],pre[v]);
	}
	if(low[u]==pre[u])
	{
		scc_cnt++;
		for(;;)
		{
			int x=s.top();s.pop();
			sccno[x]=scc_cnt;
			if(x==u)break;
		}
	}
}

void find(int n)
{
	dfs_clock=scc_cnt=0;
	memset(sccno,0,sizeof(sccno));
	memset(pre,0,sizeof(pre));
	memset(low,0,sizeof(low));
	while(!s.empty())s.pop();
	for(i=1;i<=n;i++)if(!pre[i])dfs(i);
}

void ff(int x,int y)
{
	maxx=max(maxx,y);
	for(unsigned int i=0;i<g[x].size();i++)
	{
		int v=g[x][i];
		if(!pre[v])
		{
			pre[v]=1;
			ff(v,y+tot[v]);
			pre[v]=0;
		}
	}
}	

int main()
{
	scanf("%d%d",&n,&m);
	for(i=1;i<=n;i++)scanf("%d",&a[i]);
	for(i=0;i<=n;i++)g[i].clear();
	memset(edge,0,sizeof(edge));
	for(i=1;i<=m;i++)
	{
		scanf("%d%d",&j,&k);
		g[j].push_back(k);
		edge[i][0]=j;edge[i][1]=k;
	}
	find(n);
	memset(tot,0,sizeof(tot));
	for(i=0;i<=n;i++)g[i].clear();
	for(i=1;i<=m;i++)
	{
		int p0=sccno[edge[i][0]];
		int p1=sccno[edge[i][1]];
		tot[p0]+=a[edge[i][0]];a[edge[i][0]]=0;
		tot[p1]+=a[edge[i][1]];a[edge[i][1]]=0;
		if(p0!=p1)g[p0].push_back(p1);
	}
	maxx=0;memset(pre,0,sizeof(pre));
	pre[sccno[1]]=1;ff(sccno[1],tot[sccno[1]]);
	printf("%d",maxx);
	return 0;
}
