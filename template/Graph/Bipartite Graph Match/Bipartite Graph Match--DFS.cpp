#include<bits/stdc++.h>

using namespace std;

vector<int> g[1005];
int vis[1005];
int vv[1005],res[1005];
int i,j,k,l,m,n,T;

void dfs(int x)//染色区别二分图
{
	for(unsigned int i=0;i<g[x].size();i++)
	{
		int v=g[x][i];
		if(!vis[v])
		{
			vis[v]=-vis[x];dfs(v);
			continue;
		}else if(vis[v]==vis[x])return;
	}
}

bool dfss(int x)
{
	for(unsigned int i=0;i<g[x].size();i++)
	{
		int v=g[x][i];
		if(vv[v])continue;
		vv[v]=1;//如果当前结点已被搜索过(剪枝)
		if(res[v]==0||dfss(res[v]))//寻找增广路
		{
			res[v]=x;//res表示与该点匹配的点编号
			return true;
		}
	}
	return false;
}

int main()//本程序默认二分图，非二分图会出错
{
	memset(vis,0,sizeof(vis));
	scanf("%d%d",&n,&m);
	for(i=1;i<=n;i++)g[i].clear();
	for(i=1;i<=m;i++)
	{
		scanf("%d%d",&j,&k);
		g[j].push_back(k);
		g[k].push_back(j);
	}
	for(i=1;i<=n;i++)
		if(!vis[i]){vis[i]=1;dfs(i);}//先进行黑白染色，区分开二分图
	memset(res,0,sizeof(res));k=0;
	for(i=1;i<=n;i++)//进行增广
		if(vis[i]==1)
		{
			memset(vv,0,sizeof(vv));
			if(dfss(i))k++;
		}
	printf("%d",k);//最大匹配数
	return 0;
}
