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
}a[100005];

vector<edge> q[100005];
map<string,int> my;
map<int,string> ym;

int i,j,k,l,m,n,T,o;
bool vis[100005];
int ans[100005];
int bcj[100005];

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
	for(i=1;i<=n;i++)q[i].clear();
	cin>>n;n++;
	for(i=1;i<=n;i++)bcj[i]=i;
	for(i=1;i<=n;i++)a[i].son.clear();
	my.clear();ym.clear();l=0;
	for(i=1;i<=n-1;i++)
	{
		string s1,s2;
		cin>>s1>>s2;
		j=cl(s1);k=cl(s2);
		a[j].son.push_back(k);
		vis[k]=true;
	}
	for(i=1;i<=n;i++)if(!vis[i]){o=i;break;}
	memset(vis,0,sizeof(vis));
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
	dfs(o);
	for(i=1;i<=m;i++)cout<<ym[ans[i]]<<endl;
	return 0;
}
