#include <iostream>
#include <cstring>
#include <stack>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <queue>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <deque>
#include <bitset>
#include <algorithm>
  
using namespace std;
  
const double eps=1e-10;
const double pi=3.1415926535897932384626433832795;
const double eln=2.718281828459045235360287471352;
  
#define LL long long
#define IN freopen("in.txt", "r", stdin)
#define OUT freopen("out.txt", "w", stdout)
#define scan(x) scanf("%d", &x)
#define scan2(x, y) scanf("%d%d", &x, &y)
#define scan3(x, y, z) scanf("%d%d%d", &x, &y, &z)
#define sqr(x) (x) * (x)
#define pr(x) printf("Case %d: ",x)
#define prn(x) printf("Case %d:\n",x)
#define prr(x) printf("Case #%d: ",x)
#define prrn(x) printf("Case #%d:\n",x)
#define lc o << 1
#define rc o << 1 | 1
#define lowbit(x) (x&(-x))
#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
  
const int maxn=105;
const int inf=0x3f3f3f3f;

int a[maxn];
int cost[maxn][maxn];
int lowcost[maxn],fat[maxn],maxd[maxn][maxn];
bool vis[maxn];
int i,j,k,l,m,n,T,u,v,ans,mini;

int main()
{
	scan(T);
	while(T--)
	{
		memset(lowcost,inf,sizeof(lowcost));
		memset(a,0,sizeof(a));
		memset(fat,0,sizeof(fat));
		memset(maxd,0,sizeof(maxd));
		memset(cost,inf,sizeof(cost));
		memset(vis,0,sizeof(vis));
		scan2(n,m);ans=0;
		for(i=1;i<=m;i++)
		{
			scan3(j,k,l);
			cost[j][k]=cost[k][j]=l;
		}
		vis[1]=true;a[k=1]=1;
		for(i=2;i<=n;i++){maxd[i][1]=maxd[1][i]=lowcost[i]=cost[1][i];fat[i]=1;}
		for(i=1;i<=n;i++)maxd[i][i]=cost[i][i]=0;
		for(u=1,i=1;i<=n-1;i++)
		{
			mini=inf,v=-1;
			for(j=1;j<=n;j++)
				if(!vis[j] && lowcost[j]<mini)
				{mini=lowcost[j];v=j;}
			vis[v]=true;
			ans+=mini;
			for(j=1;j<=k;j++)
				maxd[a[j]][v]=maxd[v][a[j]]=max(mini,maxd[fat[v]][a[j]]);
			a[++k]=v;
			for(j=1;j<=n;j++)
				if(!vis[j] && cost[v][j]<lowcost[j])
				{lowcost[j]=cost[v][j];fat[j]=v;}
		}
		mini=inf;
		for(i=1;i<=n-1;i++)
			for(j=i+1;j<=n;j++)
			{
				if(fat[i]==j || fat[j]==i || cost[i][j]==inf)continue;
				mini=min(mini,cost[i][j]-maxd[i][j]);
			}
		if(mini==0)printf("Not Unique!\n");else printf("%d\n",ans);
	}
    return 0;
}

