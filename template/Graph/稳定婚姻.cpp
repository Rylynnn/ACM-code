//对男性(na)最优
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

