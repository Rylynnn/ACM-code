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
i	return lowu;
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
