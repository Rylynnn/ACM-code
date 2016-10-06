#define type int//type 可选择int或者double

const type inf=2147483640;
const int maxn=1005;

int pre[maxn],id[maxn],vis[maxn];
type in[maxn];

struct edge
{
	int from,to;
	type cost;
	edge(int from=0,int to=0,type cost=0):from(from),to(to),cost(cost){}
}edg[10005];

type ZLEdmonds(int n,int m,int root)//自环在输入建图时直接忽略，如需加入，可另存
{
	type tot=0.0;
	//判断是否有树
	while(true)
	{
		for(int i=1;i<=n;i++)in[i]=inf;
		for(int i=1;i<=m;i++)
		{
			int u=edg[i].from;
			int v=edg[i].to;
			if(edg[i].cost<in[v] && u!=v){pre[v]=u;in[v]=edg[i].cost;}
		}
		for(int i=1;i<=n;i++)if(i!=root && in[i]==inf)return -1;
		//找环
		int cnt=1;
		memset(id,0,sizeof(id));
		memset(vis,0,sizeof(vis));
		in[root]=0;
		for(int i=1;i<=n;i++)//标记每个环
		{
			tot+=in[i];
			int v=i;
			while(vis[v]!=i && id[v]==0 && v!=root)
			{vis[v]=i;v=pre[v];}
			if(v!=root && id[v]==0)//缩点
			{
				for(int u=pre[v];u!=v;u=pre[u])id[u]=cnt;
				id[v]=cnt++;
			}
		}
		if(cnt==1)break;
		for(int i=1;i<=n;i++)if(id[i]==0)id[i]=cnt++;
		//建立新图
		for(int i=1;i<=m;i++)
		{
			int u=edg[i].from;
			int v=edg[i].to;
			edg[i].from=id[u];
			edg[i].to=id[v];
			if(id[u]!=id[v])edg[i].cost-=in[v];
		}
		n=cnt-1;
		root=id[root];
	}
	return tot;
}

int main()
{
	
    return 0;
}

