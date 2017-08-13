
/*****************************************
Author: lizi
Email: lzy960601@outlook.com
****************************************/
  
#include<bits/stdc++.h>
  
using namespace std;
  
const double eps=1e-10;
const double pi=3.1415926535897932384626433832795;
const double eln=2.718281828459045235360287471352;
  
#define LL long long
#define IN freopen("in.txt", "r", stdin)
#define OUT freopen("out.txt", "w", stdout)
#define scan(x) scanf("%d", &x)
#define mp make_pair
#define pb push_back
#define sqr(x) (x) * (x)
#define pr(x) printf("Case %d: ",x)
#define prn(x) printf("Case %d:\n",x)
#define prr(x) printf("Case #%d: ",x)
#define prrn(x) printf("Case #%d:\n",x)
#define lowbit(x) (x&(-x))

int st[205][10][205][10];
int a[205][205];
int lo[257];
int n,m,ans=0;

int read()
{
    char c=getchar();
    while(c<'a' || c>'z')c=getchar();
    return c-97;
}

int find(int u,int v,int x,int y)
{
    int p=lo[x-u+1],q=lo[y-v+1];
    int t=st[u][p][v][q];
    t|=st[x-(1<<p)+1][p][v][q];
    t|=st[u][p][y-(1<<q)+1][q];
    t|=st[x-(1<<p)+1][p][y-(1<<q)+1][q];
    int cnt=0;
    while(t>0)
    {
        cnt++;
        t-=lowbit(t);
        if(cnt>2)break;
    }
    return cnt;
}

int main()
{
    lo[1]=0;
    for(int i=2;i<=256;i++)lo[i]=lo[i>>1]+1;
    scanf("%d%d",&n,&m);
    for(int i=1;i<=n;i++)
        for(int j=1;j<=m;j++)
            a[i][j]=read();
    for(int i=1;i<=n;i++)
        for(int j=1;j<=m;j++)
            st[i][0][j][0]=1<<a[i][j];
    for(int mi=0;(1<<mi)<=n;mi++)
        for(int mj=0;(1<<mj)<=m;mj++)
        {
            if(mi+mj<=0)continue;
            for(int i=1;i+(1<<mi)-1<=n;i++)
                for(int j=1;j+(1<<mj)-1<=m;j++)
                    if(mi==0)st[i][mi][j][mj]|=st[i][mi][j][mj-1]|st[i][mi][j+(1<<(mj-1))][mj-1];
                    else st[i][mi][j][mj]|=st[i][mi-1][j][mj]|st[i+(1<<(mi-1))][mi-1][j][mj];
        }
    for(int i=1;i<=n;i++)
        for(int j=1;j<=m;j++)
        {
            int d=n;
            for(int k=j;k<=m;k++)
            {
                while(find(i,j,d,k)>2 && d>=i)d--;
                if(d>=i)ans+=d-i+1;
            }
        }
    printf("%d\n",ans);
    return 0;
}
