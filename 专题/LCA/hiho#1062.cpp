#include<bits/stdc++.h>

using namespace std;

struct edge
{
	int to,dist;
};

struct tree
{
	vector<int> son;
	int fat,dep;
}a[10005];

vector<edge> q[10005];
map<string,int> my;
map<int,string> ym;

int i,j,k,l,m,n,T,o;
bool vis[10005];
int ans[10005];
int bcj[10005];

int cl(string t)
{
	if(!my.count(t))
	{
		l++;
		my[t]=l;
		ym[l]=t;
		return l;
	}
	return my[t];
}

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
		if(vis[r.to]&&(!ans[r.dist]))
		{
			int zx=find(r.to);
			ans[r.dist]=zx;
		}

	}
	for(unsigned int i=0;i<a[x].son.size();i++)
	{
		int v=a[x].son[i];
		if(!vis[v])
		{
			dfs(v);
			bcj[v]=x;
		}
	}
}

int main()
{
	memset(a,sizeof(a),0);
	memset(ans,0,sizeof(ans));
	memset(vis,0,sizeof(vis));
	cin>>n;
	for(i=1;i<=500;i++)q[i].clear();
	for(i=1;i<=500;i++)bcj[i]=i;
	for(i=1;i<=500;i++)a[i].son.clear();
	my.clear();ym.clear();l=0;
	ym[10000]="-1";bcj[10000]=10000;
	for(i=1;i<=n;i++)
	{
		string s1,s2;
		cin>>s1>>s2;
		j=cl(s1);k=cl(s2);
		a[j].son.push_back(k);
		vis[k]=true;
	}
	cin>>m;
	for(i=1;i<=m;i++)
	{
		string s1,s2;
		cin>>s1>>s2;
		j=cl(s1);k=cl(s2);
		edge r;
		r.to=k;r.dist=i;
		q[j].push_back(r);
		r.to=j;r.dist=i;
		q[k].push_back(r);
	}
	memset(vis,0,sizeof(vis));
	for(i=1;i<=l;i++)if(!vis[i]&&ym.count(i))a[10000].son.push_back(i);
	dfs(10000);
	for(i=1;i<=m;i++)cout<<ym[ans[i]]<<endl;
	return 0;
}
