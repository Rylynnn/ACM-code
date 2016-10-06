#include<bits/stdc++.h>

using namespace std;

struct tree
{
	vector<int> son;
	int fat,dep;
}a[100005];

map<string,int> my;
map<int,string> ym;

int i,j,k,l,m,n,T,o,md,sta;
queue<int> q;
bool vis[100005];
int ans[100005];
int lca[100005][20];

void bfs(int x)
{
	while(!q.empty())q.pop();
	q.push(x);
	while(!q.empty())
	{
		int t=q.front();q.pop();
		a[t].dep=a[a[t].fat].dep+1;
		if(a[t].dep>md)md=a[t].dep;
		for(unsigned int i=0;i<a[t].son.size();i++)
			q.push(a[t].son[i]);
	}
}

int find(int x,int y)
{
	if(a[x].dep<a[y].dep)
	{int temp=x;x=y;y=temp;}
	while(a[x].dep>a[y].dep)
	{
		int j=0;
		while(a[lca[x][j]].dep>a[y].dep)j++;
		if(a[lca[x][j]].dep==a[y].dep){x=lca[x][j];break;}
		x=lca[x][--j];
	}
	if(x==y)return y;
	while(x!=y)
	{
		j=0;
		while(lca[x][j]!=lca[y][j])j++;
		if(j==0)break;j--;
		x=lca[x][j];y=lca[y][j];
	}
	return a[x].fat;
}

int main()
{
	memset(a,sizeof(a),0);
	memset(vis,0,sizeof(vis));
	cin>>n;
	for(i=1;i<=100005;i++)a[i].son.clear();
	my.clear();ym.clear();l=0;
	for(i=1;i<=n;i++)
	{
		string s1,s2;
		cin>>s1>>s2;
		j=cl(s1);k=cl(s2);
		a[k].fat=j;
		a[j].son.push_back(k);
		vis[k]=true;
	}
	for(i=1;i<=l;i++)if(!vis[i]){sta=i;break;}
	memset(vis,0,sizeof(vis));
	md=a[0].dep=a[0].fat=a[sta].fat=0;bfs(sta);
	for(i=1;i<=l;i++)lca[i][0]=a[i].fat;
	for(j=1;(1<<j)<=md;j++)
		for(i=1;i<=l;i++)
			lca[i][j]=lca[lca[i][j-1]][j-1];
	T=o=0;
	cin>>m;
	for(i=1;i<=m;i++)
	{
		string s1,s2;
		cin>>s1>>s2;
		j=cl(s1);k=cl(s2);
		cout<<ym[find(j,k)]<<endl;
	}
	return 0;
}
