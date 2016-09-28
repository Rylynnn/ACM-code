#include<cstdio>
#include<cstring>
#include<vector>

using namespace std;

int bcj[10005];
bool vis[10005];
int i,j,k,l,m,n,T;
vector<int> a[10005];

int find(int x)
{
	if(bcj[x]==x)return x;
	else return bcj[x]=find(bcj[x]);
}

void dfs(int x)
{
	vis[x]=true;
	if(vis[m]&&x==l)k=find(m);
	if(vis[l]&&x==m)k=find(l);
	for(unsigned int i=0;i<a[x].size();i++)
	{
		int v=a[x][i];
		if(!vis[v])
		{
			dfs(v);
			bcj[v]=x;
		}
	}
}

int main()
{
	scanf("%d",&T);
	while(T--)
	{
		memset(vis,0,sizeof(vis));
		scanf("%d",&n);
		for(i=1;i<=n;i++)bcj[i]=i;
		for(i=1;i<=n;i++)a[i].clear();
		for(i=1;i<=n-1;i++)
		{
			scanf("%d%d",&j,&k);
			a[j].push_back(k);
			vis[k]=true;
		}
		for(i=1;i<=n;i++)if(!vis[i])j=i;
		memset(vis,0,sizeof(vis));
		scanf("%d%d",&l,&m);
		dfs(j);
		printf("%d\n",k);
	}
	return 0;
}
