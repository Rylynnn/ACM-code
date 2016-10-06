#include<bits/stdc++.h>

using namespace std;

struct tree
{
	vector<int> son;
	int fat,dep;
}a[100005];

int i,j,k,l,m,n,T,o,md,sta;
bool vis[100005];
int ans[100005];
int st[100005][20];
int se[2100000],no[100005],fir[1000005];

void dfs(int x)
{
	no[++md]=x;
	int t=md;
	se[++T]=t;
	fir[x]=T;
	for(unsigned int i=0;i<a[x].son.size();i++)
	{
		dfs(a[x].son[i]);
		se[++T]=t;
	}
}

void getST()
{
	for(int i=1;i<=T;i++)st[i][0]=se[i];
	for(int j=1;(1<<j)<=T;j++)
		for(int i=1;i<=T-(1<<j)+1;i++)
			st[i][j]=min(st[i][j-1],st[i+(1<<(j-1))][j-1]);
}

int find(int x,int y)
{
	if(x>y){int temp=x;x=y;y=temp;}
	int j=0;
	while((1<<j)<=(y-x+1))j++;
	j--;
	return min(st[x][j],st[y-(1<<j)+1][j]);
}

int main()
{
	memset(a,sizeof(a),0);
	memset(vis,0,sizeof(vis));
	cin>>n;md=T=0;
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
	memset(st,0,sizeof(st));
	dfs(sta);cin>>m;
	getST();
	for(i=1;i<=m;i++)
	{
		string s1,s2;
		cin>>s1>>s2;
		j=cl(s1);k=cl(s2);
		cout<<ym[no[find(fir[j],fir[k])]]<<endl;
	}
	return 0;
}
