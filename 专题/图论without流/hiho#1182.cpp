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

const int maxn=70000;
vector<int> g[maxn];
int a[maxn];
bool b[maxn][2];
int i,j,k,l,m,n;

void dfs(int x)
{
	for(unsigned int i=0;i<g[x].size();i++)
	{
		if(b[x][i])continue;
		b[x][i]=true;
		dfs(g[x][i]);
	}
	a[++k]=x;
}

int main()
{
	scanf("%d",&n);
	m=0;
	for(i=0;i<=(1<<(n-1))-1;i++)
	{
		int p1=i<<1;
		int p2=p1^1;
		p1-=((p1>>(n-1))<<(n-1));
		p2-=((p2>>(n-1))<<(n-1));
		g[i].push_back(p1);
		g[i].push_back(p2);
	}
	k=0;dfs(0);
	if(n==1)printf("01");else for(i=1;i<=k-1;i++)printf("%d",a[i]&1);
    return 0;
}
