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
#define lowbit(x) (x&(-x))
  
const int maxn=2005;

int na[maxn][maxn],nv[maxn][maxn];
queue<int> q;
int i,j,k,m,n,T;

int main()
{
    scanf("%d",&T);
    while(T--)
    {
        scanf("%d",&n);
        memset(na,0,sizeof(na));
        memset(nv,0,sizeof(nv));
        for(i=1;i<=n;i++)
            for(j=1;j<=n;j++)scanf("%d",&na[i][j]);
        for(i=1;i<=n;i++)
            for(j=1;j<=n;j++)
            {
                scanf("%d",&m);
                nv[i][m]=j;
            }
        while(!q.empty())q.pop();
        for(i=1;i<=n;i++)q.push(i);
        while(!q.empty())
        {
            m=q.front();
            q.pop();
            k=na[m][++na[m][0]];
            if(nv[k][0]==0)
            {
                nv[k][0]=m;
                continue;
            }else
            {
                j=nv[k][0];
                if(nv[k][m]<nv[k][j])
                {
                    q.push(j);
                    nv[k][0]=m;
                    continue;
                }else q.push(m);
            }
        }
        for(i=1;i<=n;i++)na[nv[i][0]][0]=i;
        for(i=1;i<=n;i++)printf("%d\n",na[i][0]);
        if(T>0)printf("\n");
    }
    return 0;
}

