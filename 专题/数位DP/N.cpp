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

LL find(LL x)
{
    if(x==0)return 1;
    if(x<0)return 0;
    LL t,y,sum;
    y=x;t=0;sum=x/10;
    while(y>0)
    {
        t+=y%10;
        y/=10;
    }
    t-=x%10;
    for(LL i=0;i<=x%10;i++)
        if((t+i)%10==0)return sum+1;
    return sum;
}

int main()
{
    int t,T;
    scanf("%d",&T);
    for(t=1;t<=T;t++)
    {
        prr(t);
        LL x,y;
        scanf("%lld %lld",&x,&y);
        printf("%lld\n",find(y)-find(x-1));
    }
    /*while(1)
    {
        LL n;
        scanf("%lld",&n);
        printf("%lld\n",find(n));
    }*/
    return 0;
}

