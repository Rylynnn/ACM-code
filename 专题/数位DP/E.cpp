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
  
LL a[25][10];
LL b[25],c[25];
LL i,j,k,l,m,n,T;

int main()
{
    memset(a,0,sizeof(a));
    a[2][4]=1;a[2][0]=1;
    for(i=1;i<=9;i++)a[3][i]=1;
    a[3][4]=11;l=100;a[3][0]=19;
    for(i=4;i<=19;i++)
    {
        l*=10;a[i][0]=0;
        LL sum=0;
        for(j=1;j<=i-1;j++)sum+=a[j][0];
        for(j=1;j<=9;j++)a[i][j]=sum;
        a[i][4]-=a[i-1][9];
        a[i][4]+=l/10;
        for(j=1;j<=9;j++)a[i][0]+=a[i][j];
    }
    c[1]=a[1][0];
    for(i=2;i<=19;i++)c[i]=c[i-1]+a[i][0];
    scanf("%lld",&T);
    while(T--)
    {
        scanf("%lld",&n);
        m=n;k=0;l=1;
        while(m>0)
        {
            b[++k]=m%10;
            m/=10;l*=10;
        }
        LL sum=0;
        for(i=1;i<k;i++)sum+=a[i][0];
        m=n;b[k+1]=0;
        for(i=k;i>=1;i--)
        {
            m=m%l;l/=10;
            if(b[i]>1)for(j=1;j<=b[i]-1;j++)sum+=a[i][j];
            if(i<k && b[i]>0)sum+=c[i-1];
            if(b[i]==9 && b[i+1]==4)
            {
                sum+=m%l;
                sum++;
                break;
            }
        }
        printf("%lld\n",sum);
    }
    return 0;
}

