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
