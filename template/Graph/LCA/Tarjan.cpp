#include<bits/stdc++.h>

using namespace std;

struct edge
{
	int to,dist;
};

vector<edge> g[40005];
vector<edge> q[40005];

int i,j,k,l,m,n,T,o;
bool vis[40005];
int dis[40005];
int ans[205];
int bcj[40005];

int find(int x)
{
	if(bcj[x]==x)return x;
	else return bcj[x]=find(bcj[x]);
}

void dfs(int x)
{
	vis[x]=true;
	for(unsigned int i=0;i<q[x].size();i++)
	{
		edge r=q[x][i];
		if(vis[r.to])
		{
			int zx=find(r.to);
			ans[r.dist]=dis[x]+dis[r.to]-2*dis[zx];
		}
	}
	for(unsigned int i=0;i<g[x].size();i++)
	{
		edge v=g[x][i];
		if(!vis[v.to])
		{
			dis[v.to]=dis[x]+v.dist;
			dfs(v.to);
			bcj[v.to]=x;
		}
	}
}

int main()
{
	scanf("%d",&T);
	while(T--)
	{
		memset(ans,0,sizeof(ans));
		memset(vis,0,sizeof(vis));
		memset(dis,0,sizeof(dis));

		scanf("%d%d",&n,&m);
		for(i=1;i<=n;i++)q[i].clear();
		for(i=1;i<=n;i++)bcj[i]=i;
		for(i=1;i<=n;i++)g[i].clear();

		for(i=1;i<=n-1;i++)
		{
			scanf("%d%d%d",&j,&k,&l);
			edge r;
			r.to=k;r.dist=l;
			g[j].push_back(r);
			r.to=j;r.dist=l;
			g[k].push_back(r);
		}
		
		for(i=1;i<=m;i++)
		{
			scanf("%d%d",&j,&k);
			edge r;
			r.to=k;r.dist=i;
			q[j].push_back(r);
			r.to=j;r.dist=i;
			q[k].push_back(r);
		}
		memset(vis,0,sizeof(vis));
		dfs(1);
		for(i=1;i<=m;i++)printf("%d\n",ans[i]);
	}
	return 0;
}
