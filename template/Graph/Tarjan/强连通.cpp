//ÓÐÏòÍ¼

int pre[maxn],low[maxn],a[maxn],sccno[maxn],tot[maxn];
int edge[100005][2];
int dfs_clock,scc_cnt,maxx;
vector<int> g[maxn];
stack<int> s;

int n,m,p0,p1,i,j,k,l;

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
