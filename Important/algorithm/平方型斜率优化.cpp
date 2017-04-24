const int maxn=500005;
LL dp[maxn];
LL a[maxn];
int q[maxn];
int n,i,j;
int h,t;
LL m;

LL ga(int x,int y)
{
    LL t=a[y]-a[x];
    t-=m;
    return dp[x]+t*t;
}

LL gu(int x,int y)
{
    LL t=(a[y]-a[x])*m;
    t*=2;
    return dp[y]-dp[x]+a[y]*a[y]-a[x]*a[x]+t;
}

LL gd(int x,int y)
{
    return 2*(a[y]-a[x]);
}

int main()
{
    while(scanf("%d%lld",&n,&m)==2)
    {
        m++;
        for(i=1;i<=n;i++)scanf("%lld",&a[i]),a[i]++;
        for(a[0]=0,i=1;i<=n;i++)a[i]+=a[i-1];
        dp[0]=0;
        h=t=0;q[0]=0;
        for(i=1;i<=n;i++)
        {
            while(t-h>0 && gu(q[h],q[h+1])<=a[i]*gd(q[h],q[h+1]))h++;
            dp[i]=ga(q[h],i);
            while(t-h>0 && gu(q[t-1],q[t])*gd(q[t],i)>=gu(q[t],i)*gd(q[t-1],q[t]))t--;
            q[++t]=i;
        }
        printf("%lld\n",dp[n]);
    }
    return 0;
}

