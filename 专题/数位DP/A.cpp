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
  
char s[6];
int a[6];
int c[30][30];
int ans,i,j,k,l,m,n;

int main()
{
    memset(c,0,sizeof(c));
    c[1][0]=c[1][1]=1;
    for(i=2;i<=26;i++)
    {
        c[i][0]=c[i][i]=1;
        for(j=1;j<=min(5,i-1);j++)
            c[i][j]=c[i-1][j]+c[i-1][j-1];
    }
    while(gets(s)!=NULL)
    {
        bool f=true;
        ans=0;
        n=strlen(s);
        for(i=0;i<n-1;i++)
            if(s[i+1]<=s[i]){f=false;break;}
        if(!f){printf("0\n");continue;}
        for(i=1;i<=n-1;i++)ans+=c[26][i];
        for(i=1;i<=n;i++)a[i]=s[i-1]-96;
        a[0]=0;
        for(i=1;i<=n;i++)
            if(a[i]-a[i-1]>1)
            {
                for(j=a[i-1]+1;j<a[i];j++)
                    ans+=c[26-j][n-i];
            }
        printf("%d\n",ans+1);
    }
    return 0;
}

