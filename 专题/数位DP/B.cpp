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
  
int a[10];
int s[1000001];
int i,j,k,l,m,n;

int main()
{
    l=0;s[0]=0;
    for(i=1;i<=1000000;i++)
    {
        s[i]=s[i-1];
        j=i;k=0;
        bool f=true;
        while(j>0)
        {
            a[++k]=j%10;
            j/=10;
        }
        for(j=1;j<=k;j++)if(a[j]==4){f=false;break;}
        if(!f)continue;
        for(j=1;j<=k-1;j++)if(a[j]==2 && a[j+1]==6){f=false;break;}
        if(f)s[i]++;
    }
    while(scan2(n,m)==2 && n+m>0)
    {
            printf("%d\n",s[m]-s[n-1]);
    }
    return 0;
}

