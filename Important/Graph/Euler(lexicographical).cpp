const int maxn=10005;
vector<int> g[maxn];
int du[maxn];
bool mp[maxn][maxn];
int ans[10*maxn];
int i,j,k,l,m,n,t;

void dfs(int x)
{
    for(unsigned int i=0;i<g[x].size();i++)
    {
        int v=g[x][i];
        if(!mp[x][v])
        {
            mp[x][v]=mp[v][x]=true;
            dfs(v);
        }
    }
    ans[++t]=x;
}

struct init
{
    int x,y;
    bool operator < (const struct init p)const
    {
        if(x<p.x)return true;
        if(x>p.x)return false;
        return y<p.y;
    }
}dat[100005];

bool check()
{
    bool vis[maxn];
    memset(vis,0,sizeof(vis));
    queue<int> q;
    while(!q.empty())q.pop();
    q.push(1);vis[1]=true;
    while(!q.empty())
    {
        int now=q.front();q.pop();
        for(unsigned int i=0;i<g[now].size();i++)
        {
            if(vis[g[now][i]])continue;
            q.push(g[now][i]);
            vis[g[now][i]]=true;
        }
    }
    for(int i=1;i<=n;i++)
        if(!vis[i])return false;
    return true;
}

int main()
{
	//未判断图是否连通，为严谨可加入bfs检查
	scan2(n,m);
	for(i=1;i<=n;i++)g[i].clear();
	memset(du,0,sizeof(du));
    for(i=1;i<=m;i++)
    {
        scanf("%d%d",&dat[i].x,&dat[i].y);
        if(dat[i].x>dat[i].y)swap(dat[i].x,dat[i].y);
    }
    sort(dat+1,dat+m+1);
	for(i=1;i<=m;i++)
	{
	    j=dat[i].x;k=dat[i].y;	
		du[j]++;du[k]++;
		g[j].push_back(k);
		g[k].push_back(j);
	}
	int tot=0,st=0;
	for(i=n;i>=1;i--)
	{
		if(du[i]==0){tot=3;break;}
		if(du[i]%2==1){tot++;st=i;}
	}
	if(st==0)st=1;
	if(tot>2 || !check())printf("-1");else dfs(st);
    for(i=t;i>=1;i--)
    {
        printf("%d",ans[i]);
        if(i>1)printf(" ");
    }
    return 0;
}
