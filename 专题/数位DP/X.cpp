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
  
int a[40];
LL c[40][40];
int n,m;

void calc()
{
    memset(c,0,sizeof(c));
    c[0][0]=1;
    c[1][0]=c[1][1]=1;
    for(int i=2;i<=30;i++)
    {
        c[i][0]=c[i][i]=1;
        for(int j=1;j<=i-1;j++)c[i][j]=c[i-1][j]+c[i-1][j-1];
    }
}

LL find(int x)
{
    if(x==0)return 1;
    LL y=x;
    int k=0,j,t;
    while(y>0)
    {
        a[++k]=y%2;
        y/=2;
    }
    int odd,even;
    even=0;
    odd=1;
    LL sum=1;
    for(int i=k-1;i>=1;i--)
    {
        t=i-1;j=t;
        while(j>=t-j+1 && j>=0)
        {
            sum+=c[t][j];
            j--;
        }
    }
    for(int i=k-1;i>=1;i--)
    {
        if(a[i]==0){even++;continue;}
        t=i-1;j=t;
        while(even+j+1>=odd+(t-j) && j>=0)
        {
            sum+=c[t][j];
            j--;
        }
        odd++;
    }
    if(even>=odd)sum++;
    return sum;
}

int main()
{
    calc();
    scanf("%d %d",&n,&m);
    printf("%lld",find(m)-find(n-1));
    return 0;
}
