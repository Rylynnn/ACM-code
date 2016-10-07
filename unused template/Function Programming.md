# Graph

***

###Dijkstra(nlogn)

***

```cpp
#define pii pair<int, int>
priority_queue<pii, vector<pii>, greater<pii> > heap; 
struct edge
{
	int v, l, next;
}e[13007];
int n, m, S, T;
int tot=2, dist[2502], head[2502];
bool visited[2502];
void addedge(int x,int y,int z)
{
	e[tot].v=y, e[tot].l=z, e[tot].next=head[x], head[x]=tot++;
	e[tot].v=x, e[tot].l=z, e[tot].next=head[y], head[y]=tot++;
}
void dij(int x)
{
	if (x == T) return ;
	visited[x] = true;
	for (int p = head[x]; p; p = e[p].next)
		if (!visited[e[p].v] && dist[e[p].v] > (dist[x] + e[p].l))
			dist[e[p].v] = dist[x] + e[p].l, heap.push(make_pair(dist[e[p].v], e[p].v));
	while (!heap.empty() && visited[heap.top().second])
		heap.pop();
	dij(heap.top().second);
}		
int main()
{
	scanf("%d%d%d%d",&n,&m,&S,&T);
	int x, y, z;
	for (int i=1;i<=m;i++)
		scanf("%d%d%d",&x,&y,&z), addedge(x, y, z);
	for (int i=1;i<=n;i++)
		dist[i]=0x7ffffff;
	dist[S]=0;
	dij(S);
	printf("%d\n",dist[T]);
	return 0;
}
```

###K-th Shortest Path

```cpp
//若需要记录路径，在node里单独加一个变量，并在A*算法中记录当前节点(k)的第cnt[k]的父亲即可
struct node
{
         int v,c;
         node(int v,int c):v(v),c(c){}
         inline bool operator<(const node &b) const//用于优先队列先出的条件
         {
                   return c+dis[v]>b.c+dis[b.v];
         }
};
vector<node> map1[MAXN];//用于dijkstra算法
vector<node> map2[MAXN];//用于A_star算法
void dijkstra()
{
         int i,find[MAXN],v;
         for(i=1;i<=n;i++)dis[i]=INF;
         memset(find,0,sizeof(find));
         priority_queue<node> heap;
         dis[t]=0;
         heap.push(node(t,0));
         while(!heap.empty())
         {
                   v=heap.top().v;
                   heap.pop();
                   if(find[v])continue;
                   find[v]=1;
                   for(i=0;i<map1[v].size();i++)
                            if(!find[map1[v][i].v] && dis[v]+map1[v][i].c<dis[map1[v][i].v])
                            {
                                     dis[map1[v][i].v]=dis[v]+map1[v][i].c;
                                     heap.push(node(map1[v][i].v,dis[map1[v][i].v]));
                            }
         }
}
int A_star()
{
         int i,cnt[MAXN],v,g;
         if(dis[s]==INF)return -1;
         priority_queue<node> heap;
         memset(cnt,0,sizeof(cnt));
         heap.push(node(s,0));//0是g（x）
         while(!heap.empty())
         {
                   v=heap.top().v;
                   g=heap.top().c;
                   heap.pop();
                   cnt[v]++;
                   if(cnt[t]==k)return g;
                   if(cnt[v]>k)continue;
                   for(i=0;i<map2[v].size();i++)
                            heap.push(node(map2[v][i].v,g+map2[v][i].c));
         }
         return -1;
}
int main()
{
         int i,u,v,c;
         cin>>n>>m;
         for(i=0;i<m;i++)
         {
                   cin>>u>>v>>c;
                   map2[u].push_back(node(v,c));
                   map1[v].push_back(node(u,c));//反向储存求各节点到目标节点的最短距离
         }
         cin>>s>>t>>k;
         if(s==t)k++;
         dijkstra();
         int ans=A_star();
         cout<<ans<<endl;
         return 0;
}
```

###Hamilton

```cpp
int map[Max][Max];
int ans[Max];
bool vis[Max];

//ans数组的index 
int index;
int n, m;
int s, t;

void init()
{
    for (int i = 0; i < Max; ++i)
        for (int j = 0; j < Max; ++j)
            if (i == j)
                map[i][j] = 0;
            else
                map[i][j] = 1;
    
    memset(ans, 0, sizeof(ans));
    memset(vis, 0 , sizeof(vis));
    index = 0;
}

void reverse(int a, int b)
{
    while (a < b)
    {
        swap(ans[a], ans[b]);
        a++; 
        b--;
    }
}

void expand()
{
    while (true)
    {
        int i;
        for (i = 1; i <= n; ++i) 
        {
            if (!vis[i] && map[i][t])//未被访问且与t相连 
            {
                ans[index++] = i;
                vis[i] = true;
                t = i;       
                break;
            }
        }
        if (i > n) break;//无法扩展 
    }
}

void Hamilton()
{
    //初始化s = 1
    s = 1; 
    
    //取任意连接s的点 
    for (int i = 1; i <= n; ++i)
    {
        if (map[i][s])
        {
            t = i;       
            break;
        }
    }
    vis[s] = true;
    vis[t] = true;
    ans[index++] = s;
    ans[index++] = t;
    
    while (true)
    {
        //从t向外扩展 
        expand(); 
        
        //t扩展完毕，倒置ans并交换s,t 
        reverse(0, index-1);
        
        swap(s, t);
        
        //从另一头，t(原来的s)继续扩展 
        expand();
        
        //若s,t不相连,处理成相连 
        if (!map[s][t])
        {
            //在ans[1]到ans[index-2]中寻找两个相邻的且与st同时相连的点（必存在） 因为涉及i+1所以i < index-2 
            for (int i = 1; i < index-2; ++i) 
            {
                if (map[ans[i+1]][s] && map[ans[i]][t])
                {
                    reverse(i+1, index-1);//倒置ans[i+1]到ans[index-1] 
                    t = ans[index-1];//更新t 
                    break;
                }
            }
        }
        
        //若ans元素有n个，说明算法完成 
        if (index == n) return;
        
        //若ans元素不满n个，ans[]中寻找与未被遍历过的点相连的点，但这一点必定不是s,t.因为s,t已经遍历到无法遍历才能走到这一步 
        for (int j = 1; j <= n; ++j)
        {
            if (!vis[j])
            {
                int i;
                for (i = 1; i < index-1; ++i)//排除st 
                {
                    if (map[ans[i]][j])
                    {
                        s = ans[i-1];
                        t = j;
                        reverse(0, i-1);
                        reverse(i,index-1);
                        ans[index++] = j;
                        vis[j] = true;       
                        break;
                    }
                }       
                if (map[ans[i]][j])break;//记得有2个循环，要break两次 
            }
        }
        //继续返回，从t扩展。。 
    }
}
```

###Euler(Fleury)

```cpp
const int maxn=10005;
int stac[maxn],sta;//栈深度最大可能为边数+1
struct edge
{
	int p,q;
	edge(int p=0,int q=0):p(p),q(q){}
}edg[2*maxn];
bool used[2*maxn];
vector<int> g[maxn];
int du[maxn];
int i,j,k,l,m,n;

int other(int num,int x)
{
	return x==edg[num].p?edg[num].q:edg[num].p;
}

void dfs(int x)
{
	stac[++sta]=x;
	for(unsigned int i=0;i<g[x].size();i++)
	{
		if(used[g[x][i]])continue;
		used[g[x][i]]=true;
		dfs(other(g[x][i],x));
		break;
	}
}

void Fleury(int x)
{
	sta=1;stac[sta]=x;
	while(sta>=1)
	{
		x=stac[sta];
		bool f=false;
		for(unsigned int i=0;i<g[x].size();i++)
		{
			if(!used[g[x][i]]){f=true;break;}
		}
		if(!f)printf("%d ",stac[sta--]);
		else 
		{
			sta--;
			dfs(stac[sta+1]);
		}
	}
}

int main()
{
	//未判断图是否连通，为严谨可加入bfs检查
	scan2(n,m);
	memset(edg,0,sizeof(edg));
	for(i=1;i<=n;i++)g[i].clear();
	memset(used,0,sizeof(used));
	memset(du,0,sizeof(du));
	for(i=1;i<=m;i++)
	{
		scan2(j,k);
		edg[i]=edge(j,k);
		du[j]++;du[k]++;
		g[j].push_back(i);
		g[k].push_back(i);
	}
	int tot=0,st=0;
	for(i=1;i<=n;i++)
	{
		if(du[i]==0){tot=3;break;}
		if(du[i]%2==1){tot++;st=i;}
	}
	if(st==0)st=1;
	if(tot>2)printf("HeHeDa!");else Fleury(st);
    return 0;
}
```

***

## Tarjan

```cpp
void tarjan(int n)
{
	dfs_clock=bcc_cnt=0;
	memset(pre,0,sizeof(pre));
	memset(bcc,0,sizeof(bcc));
	memset(p,0x3f3f3f3f,sizeof(p));//点-双
	s[0]=0;//Stack
	for(int i=1;i<=n;i++)if(!pre[i])dfs(i,-1);
}
```

###割点+桥

```cpp
int dfs(int u,int fat)
{
	int lowu=pre[u]=++dfs_clock;
	for(unsigned int i=0;i<g[u].size();i++)
	{
		int side=g[u][i];
		int v=other(side,u);
		if(!pre[v])
		{
			int lowv=dfs(v,u);
			lowu=min(lowu,lowv);
			if(lowv>pre[u])bri[side]=true;
		}else if(v!=fat)lowu=min(lowu,pre[v]);
	}
	return lowu;
}
```

###点-双连通

```cpp
int dfs(int u,int fa)
{
	int lowu=pre[u]=++dfs_clock;
	for(unsigned int i=0;i<g[u].size();i++)
	{
		int side=g[u][i];
		int v=other(side,u);
		if(!pre[v])
		{
			s[++s[0]]=side;
			int lowv=dfs(v,u);
			lowu=min(lowu,lowv);
			if(lowv>=pre[u])
			{
				bcc_cnt++;
				for(;;)
				{
					int x=s[s[0]];s[0]--;
					bcc[x]=bcc_cnt;
					p[bcc_cnt]=min(p[bcc_cnt],x);
					if(x==side)break;
				}
			}
		}else if(pre[v]<pre[u] && v!=fa)
		{
			s[++s[0]]=side;lowu=min(lowu,pre[v]);
		}
	}
	return lowu;
}
```

###边-双连通

```cpp
//DFS不走桥即可
```

###SCC(强连通分量)

```cpp
void dfs(int u)
{
	pre[u]=low[u]=++dfs_clock;
	s.push(u);
	for(unsigned int i=0;i<g[u].size();i++)
	{
		int v=g[u][i];
		if(!pre[v])
		{
			dfs(v);
			low[u]=min(low[u],low[v]);
		}else if(!sccno[v])low[u]=min(low[u],pre[v]);
	}
	if(low[u]==pre[u])
	{
		scc_cnt++;
		for(;;)
		{
			int x=s.top();s.pop();
			sccno[x]=scc_cnt;
			if(x==u)break;
		}
	}
}
```

***

###次小生成树(Using Prim)

```cpp
for(i=1;i<=m;i++)
{
	scan3(j,k,l);
	cost[j][k]=cost[k][j]=l;
}
vis[1]=true;a[k=1]=1;
for(i=2;i<=n;i++){maxd[i][1]=maxd[1][i]=lowcost[i]=cost[1][i];fat[i]=1;}
for(i=1;i<=n;i++)maxd[i][i]=cost[i][i]=0;
for(u=1,i=1;i<=n-1;i++)
{
	mini=inf,v=-1;
	for(j=1;j<=n;j++)
		if(!vis[j] && lowcost[j]<mini)
		{mini=lowcost[j];v=j;}
	vis[v]=true;
	ans+=mini;
	for(j=1;j<=k;j++)
		maxd[a[j]][v]=maxd[v][a[j]]=max(mini,maxd[fat[v]][a[j]]);
	a[++k]=v;
	for(j=1;j<=n;j++)
		if(!vis[j] && cost[v][j]<lowcost[j])
		{lowcost[j]=cost[v][j];fat[j]=v;}
}
mini=inf;
for(i=1;i<=n-1;i++)
	for(j=i+1;j<=n;j++)
	{
		if(fat[i]==j || fat[j]==i || cost[i][j]==inf)continue;
		mini=min(mini,cost[i][j]-maxd[i][j]);
	}
if(mini==0)printf("Not Unique!\n");else printf("%d\n",ans);
//次小值为ans+mini
```

###最小树形图

```cpp
type ZLEdmonds(int n,int m,int root)//自环在输入建图时直接忽略，如需加入，可另存
{
	type tot=0.0;
	//判断是否有树
	while(true)
	{
		for(int i=1;i<=n;i++)in[i]=inf;
		for(int i=1;i<=m;i++)
		{
			int u=edg[i].from;
			int v=edg[i].to;
			if(edg[i].cost<in[v] && u!=v){pre[v]=u;in[v]=edg[i].cost;}
		}
		for(int i=1;i<=n;i++)if(i!=root && in[i]==inf)return -1;
		//找环
		int cnt=1;
		memset(id,0,sizeof(id));
		memset(vis,0,sizeof(vis));
		in[root]=0;
		for(int i=1;i<=n;i++)//标记每个环
		{
			tot+=in[i];
			int v=i;
			while(vis[v]!=i && id[v]==0 && v!=root)
			{vis[v]=i;v=pre[v];}
			if(v!=root && id[v]==0)//缩点
			{
				for(int u=pre[v];u!=v;u=pre[u])id[u]=cnt;
				id[v]=cnt++;
			}
		}
		if(cnt==1)break;
		for(int i=1;i<=n;i++)if(id[i]==0)id[i]=cnt++;
		//建立新图
		for(int i=1;i<=m;i++)
		{
			int u=edg[i].from;
			int v=edg[i].to;
			edg[i].from=id[u];
			edg[i].to=id[v];
			if(id[u]!=id[v])edg[i].cost-=in[v];
		}
		n=cnt-1;
		root=id[root];
	}
	return tot;
}//type=int/double/long long int
```

###0-1分数规划

```cpp
struct hei
{
	double num;
	int pos;
	hei(double num=0,int pos=0):num(num),pos(pos){}
	bool operator < (struct hei p)const
	{return num>p.num;}
}d[maxn];
double p,q,ans,l;
int i,n,m;

int main()
{
	while(scan2(n,m)==2 && n+m>0)
	{
		m=n-m;
		for(i=1;i<=n;i++)scanf("%lf",&a[i]);
		for(i=1;i<=n;i++)scanf("%lf",&b[i]);
		l=0;
		while(true)
		{
			ans=l;
			for(i=1;i<=n;i++)d[i]=hei(a[i]-ans*b[i],i);
			sort(d+1,d+n+1);
			double fz,fm;
			fz=fm=0.0;
			for(i=1;i<=m;i++)
			{
				fz+=a[d[i].pos];
				fm+=b[d[i].pos];
			}
			l=fz/fm;
			if(fabs(ans-l)<eps)break;
		}
		printf("%.0f\n",100.0*ans);
	}
    return 0;
}
```

###稳定婚姻问题(Best For Na)

```cpp
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
```

***

## LCA

###Tarjan

```cpp
int find(int x)
{
	if(bcj[x]==x)return x;
	else return bcj[x]=find(bcj[x]);
}

void dfs(int x)
{
	vis[x]=true;
	for(unsigned int i=0;i<q[x].size();i++)
	{
		edge r=q[x][i];
		if(vis[r.to]&&(!ans[r.dist]))
		{
			int zx=find(r.to);
			ans[r.dist]=zx;
		}

	}
	for(unsigned int i=0;i<a[x].son.size();i++)
	{
		int v=a[x].son[i];
		if(!vis[v])
		{
			dfs(v);
			bcj[v]=x;
		}
	}
}
```

###ST

```cpp
void dfs(int x)
{
	no[++md]=x;
	int t=md;
	se[++T]=t;
	fir[x]=T;
	for(unsigned int i=0;i<a[x].son.size();i++)
	{
		dfs(a[x].son[i]);
		se[++T]=t;
	}
}

void getST()
{
	for(int i=1;i<=T;i++)st[i][0]=se[i];
	for(int j=1;(1<<j)<=T;j++)
		for(int i=1;i<=T-(1<<j)+1;i++)
			st[i][j]=min(st[i][j-1],st[i+(1<<(j-1))][j-1]);
}

int find(int x,int y)
{
	if(x>y){int temp=x;x=y;y=temp;}
	int j=0;
	while((1<<j)<=(y-x+1))j++;
	j--;
	return min(st[x][j],st[y-(1<<j)+1][j]);
}
```

###倍增

```cpp
void bfs(int x)
{
	while(!q.empty())q.pop();
	q.push(x);
	while(!q.empty())
	{
		int t=q.front();q.pop();
		a[t].dep=a[a[t].fat].dep+1;
		if(a[t].dep>md)md=a[t].dep;
		for(unsigned int i=0;i<a[t].son.size();i++)
			q.push(a[t].son[i]);
	}
}

int find(int x,int y)
{
	if(a[x].dep<a[y].dep)
	{int temp=x;x=y;y=temp;}
	while(a[x].dep>a[y].dep)
	{
		int j=0;
		while(a[lca[x][j]].dep>a[y].dep)j++;
		if(a[lca[x][j]].dep==a[y].dep){x=lca[x][j];break;}
		x=lca[x][--j];
	}
	if(x==y)return y;
	while(x!=y)
	{
		j=0;
		while(lca[x][j]!=lca[y][j])j++;
		if(j==0)break;j--;
		x=lca[x][j];y=lca[y][j];
	}
	return a[x].fat;
}

//Part of Main
md=a[0].dep=a[0].fat=a[sta].fat=0;bfs(sta);
for(i=1;i<=l;i++)lca[i][0]=a[i].fat;
for(j=1;(1<<j)<=md;j++)
	for(i=1;i<=l;i++)
		lca[i][j]=lca[lca[i][j-1]][j-1];
```

***

##Flow

###Dinic

```cpp
struct edge
{
	int s,t,v,c;
	edge* next;
}mem[maxm],*head[maxn],*prev[maxn];
queue<int> q;
int cnt=-1,n;
int dis[maxn];
int S,T;
void add_edge(int s,int t,int v,int c)
{
	mem[++cnt].s=s;mem[cnt].t=t;mem[cnt].v=v;mem[cnt].c=c;mem[cnt].next=head[s];head[s]=mem+cnt;
	mem[++cnt].s=t;mem[cnt].t=s;mem[cnt].v=0;mem[cnt].c=-c;mem[cnt].next=head[t];head[t]=mem+cnt;
}
bool bfs()
{
	for (int i=0;i<=n;i++) dis[i]=INF;
	q.push(S);dis[S]=0;
	while(!q.empty())
	{
		for (edge *it=head[q.front()];it;it=it->next)
		if (it->v&&dis[q.front()]+it->c<dis[it->t])
		{
			dis[it->t]=dis[q.front()]+it->c;
			prev[it->t]=it;
			q.push(it->t);
		}
		q.pop();
	}
	return (dis[T]!=INF);
}
int cost=0;
int dinic()
{
	int flow=0;
	while(bfs())
	{
		int augflow=INF,tmpcost=0;
		for (edge* it=prev[T];it;it=prev[it->s])
		{
			augflow=min(augflow,it->v);
			tmpcost+=it->c;
		}
		for (edge* it=prev[T];it;it=prev[it->s])
		{
			it->v-=augflow;
			rever(it)->v+=augflow;
		}
		flow+=augflow;cost+=augflow*tmpcost;
	}
	return flow;
}
```

###ISAP

```cpp
struct edge
{
	int s,t,v;
	edge* next;
}mem[maxm*2],*head[maxn];
int cnt=-1;
void add_edge(int s,int t,int v)
{
	mem[++cnt].s=s;mem[cnt].t=t;mem[cnt].v=v;mem[cnt].next=head[s];head[s]=mem+cnt;
	mem[++cnt].s=t;mem[cnt].t=s;mem[cnt].v=0;mem[cnt].next=head[t];head[t]=mem+cnt;
}
int n,m;

int S,T;
int numbs[maxn];
int d[maxn];					
edge* cur[maxn],*revpath[maxn];	

void bfs()
{
	queue<int> q;
	while(!q.empty()) q.pop();
	for (int i=1;i<=n;i++) d[i]=maxn-1;				//由初始下标决定01
	d[T]=0;q.push(T);
	while(!q.empty())
	{
		int u=q.front();
		q.pop();
		for (edge* it=head[u];it;it=it->next)
		{
			edge *now=rever(it);
			if (now->v==0||d[now->s]<n) continue;
			d[now->s]=d[u]+1;
			q.push(now->s);
		}
	}
	memset(numbs,0,sizeof(numbs));
	for (int i=1;i<=n;i++) numbs[d[i]]++;			//由初始下标决定01			
}

int isap()
{
	int flow=0;
	for (int i=1;i<=n;i++) cur[i]=head[i];			//由初始下标决定01
	int u=S;
	while(d[S]<n)
	{
		if (u==T) 
		{
			int augflow=2147483647;
			for (int i=S;i!=T;i=cur[i]->t)
				augflow=min(augflow,cur[i]->v);
			for (int i=S;i!=T;i=cur[i]->t)
			{
				cur[i]->v-=augflow;
				rever(cur[i])->v+=augflow;
			}
			flow+=augflow;u=S;
		}
		edge *e;
		for (e=cur[u];e;e=e->next)
			if (e->v&&d[u]==(d[e->t]+1)) break;
		if (e)
		{
			cur[u]=e;
			revpath[e->t]=rever(e);
			u=e->t;
		}
		else
		{
			numbs[d[u]]--;
			if (numbs[d[u]]==0) break;
			cur[u]=head[u];
			int mindist=n;
			for (edge* it=head[u];it;it=it->next)
				if (it->v) mindist=min(mindist,d[it->t]);
			d[u]=mindist+1;
			numbs[d[u]]++;
			if (u!=S) u=revpath[u]->t;
		}
	}
	return flow;
}
```



***

***

***



# Data Structures

***

***

***

## Tree

###Cartesian_Tree

```cpp
/*
@title: Cartesian Tree 笛卡尔树
@description:
    Cartesian Tree 笛卡尔树
    可以实现线性时间内建立具有BST性质的树
@structure:
    CartesianTreeNode:
        parent:　父指针
        l: 左孩子指针
        r: 右孩子指针
@arguments:
    BuildFromArray:
        value: 源数组
        N: 数组大小
        index: 源数组的逆映射数组
        tree: 目标建树数组内存首地址
        stack: 堆栈空间
@performance:
    BuildFromArray:
        Time: O(N)
        Space: O(N)
@dependence: null
@range:
    for i in [0, N)
    value[i] in [0, N)
    index[i] in [0, N)
    |value| = |index| = |tree| = |stack| = N
@note:
    value 与 index 互为逆映射故满足双射性质
        index[value[i]] == i
        value[index[i]] == i
    index 无须在函数外初始化，建树过程可以计算 index
    stack 无须在函数外初始化，但建树过程对 stack 有污染
    最后结束迭代的时候栈底一定为 value[0]
    笛卡尔树的树根一定为 value[0]
    因此笛卡尔树的 parent 不一定要保存，仅保存孩子指针也可以完成遍历
*/

struct CartesianTreeNode {
  int parent, left, right;
};

void BuildFromArray(int *value, int N, int *index, CartesianTreeNode *tree,
                    int *stack) {
  // 计算逆映射
  for (int i = 0; i < N; i++) {
    index[value[i]] = i;
  }
  // 初始化节点
  for (int i = 0; i < N; i++) {
    tree[i].parent = tree[i].left = tree[i].right = -1;
  }
  int size = 0; // 初始化清空栈
  for (int i = 0; i < N; i++) {
    int nextSize = size;
    // 维护单调栈
    while (nextSize > 0 && index[stack[nextSize - 1]] > index[i]) {
      nextSize--;
    }
    // 下面两个 if 语句块的顺序可变
    if (nextSize > 0) { // 栈中有元素
      // 当前元素为栈顶元素的右孩子
      int top = stack[nextSize - 1];
      tree[i].parent = top;
      tree[top].right = i;
    }
    if (nextSize < size) { // 弹过栈
      // 最后出栈的元素为当前元素的左孩子
      int lastPop = stack[nextSize];
      tree[lastPop].parent = i;
      tree[i].left = lastPop;
    }
    stack[nextSize++] = i; // 入栈
    size = nextSize;       // 更新栈大小
  }
}
```

###Segment Tree

```cpp
void update(int x)
{
	if(a[x].ls+a[x].rs==0)return;
	a[x].mi=min(a[a[x].ls].mi,a[a[x].rs].mi);
}

void downdate(int x)
{
	if(a[x].ls+a[x].rs==0)return;
	if(a[x].add==0)return;
	a[a[x].ls].add+=a[x].add;
	a[a[x].ls].mi+=a[x].add;
	a[a[x].rs].add+=a[x].add;
	a[a[x].rs].mi+=a[x].add;
	a[x].add=0;
}

void mt(int l,int r)
{
	if(l==r)
	{
		a[k].ll=a[k].rr=l;
		a[k].mi=c[l];
		a[k].ls=a[k].rs=0;
		return;
	}
	int t=k;
	int mid=(l+r)>>1;
	a[t].ll=l;a[t].rr=r;
	k++;a[t].ls=k;mt(l,mid);
	k++;a[t].rs=k;mt(mid+1,r);
	update(t);
}

void add(int l,int r,int nu,int x)
{
	if(a[x].ll==l && a[x].rr==r)
	{
		a[x].add+=nu;
		a[x].mi+=nu;
		return;
	}
	downdate(x);
	int mid=(a[x].ll+a[x].rr)>>1;
	if(mid<l)add(l,r,nu,a[x].rs);
	else if(mid>=r)add(l,r,nu,a[x].ls);
	else
	{
		add(l,mid,nu,a[x].ls);
		add(mid+1,r,nu,a[x].rs);
	}
	update(x);
}

int find(int l,int r,int x)
{
	if(a[x].ll==l && a[x].rr==r)return a[x].mi;
	downdate(x);
	update(x);
	int mid=(a[x].ll+a[x].rr)>>1;
	if(mid<l)return find(l,r,a[x].rs);
	if(mid>=r)return find(l,r,a[x].ls);
	return min(find(l,mid,a[x].ls),find(mid+1,r,a[x].rs));
}
```

###Treap

```cpp
void update(int k)//更新结点信息
{
    tr[k].size=tr[tr[k].l].size+tr[tr[k].r].size+tr[k].w;
}

void rturn(int &k)
{
    int t=tr[k].l;tr[k].l=tr[t].r;tr[t].r=k;
    tr[t].size=tr[k].size;update(k);k=t;
}

void lturn(int &k)
{
    int t=tr[k].r;tr[k].r=tr[t].l;tr[t].l=k;
    tr[t].size=tr[k].size;update(k);k=t;
}

void insert(int &k,int x)
{
    if(k==0)
    {
        size++;k=size;
        tr[k].size=tr[k].w=1;tr[k].v=x;tr[k].rnd=rand();
        return;
    }
    tr[k].size++;
    if(tr[k].v==x)tr[k].w++;
    else if(x>tr[k].v)
    {
        insert(tr[k].r,x);
        if(tr[tr[k].r].rnd<tr[k].rnd)lturn(k);
    }
    else 
    {
        insert(tr[k].l,x);
        if(tr[tr[k].l].rnd<tr[k].rnd)rturn(k);
    } 
}

void del(int &k,int x)
{
    if(k==0)return; 
    if(tr[k].v==x)
    {
        if(tr[k].w>1)
        {
            tr[k].w--;tr[k].size--;return;
        }
        if(tr[k].l*tr[k].r==0)k=tr[k].l+tr[k].r;
        else if(tr[tr[k].l].rnd<tr[tr[k].r].rnd)
            rturn(k),del(k,x);
        else lturn(k),del(k,x);
    }
    else if(x>tr[k].v)
        tr[k].size--,del(tr[k].r,x);
    else tr[k].size--,del(tr[k].l,x);
}

int query_rank(int k,int x)
{
    if(k==0)return 0;
    if(tr[k].v==x)return tr[tr[k].l].size+1;
    else if(x>tr[k].v)
        return tr[tr[k].l].size+tr[k].w+query_rank(tr[k].r,x);
    else return query_rank(tr[k].l,x);
}

int query_num(int k,int x)
{
    if(k==0)return 0;
    if(x<=tr[tr[k].l].size)
        return query_num(tr[k].l,x);
    else if(x>tr[tr[k].l].size+tr[k].w)
        return query_num(tr[k].r,x-tr[tr[k].l].size-tr[k].w);
    else return tr[k].v;
}

void query_pro(int k,int x)
{
    if(k==0)return;
    if(tr[k].v<x)
    {
        ans=k;query_pro(tr[k].r,x);
    }
    else query_pro(tr[k].l,x);
}

void query_sub(int k,int x)
{
    if(k==0)return;
    if(tr[k].v>x)
    {
        ans=k;query_sub(tr[k].l,x);
    }
    else query_sub(tr[k].r,x);
}
```

###Splay

```cpp
struct tree
{
	int key,size,le,ri,add,rev,min,pre;
}a[maxn];
int n,T,node;
int s[maxn];

void pushdown(int cur)
{
	int ls=a[cur].le,rs=a[cur].ri;
	if(a[cur].add>0)
	{
		a[ls].add+=a[cur].add;
		a[rs].add+=a[cur].add;
		a[ls].key+=a[cur].add;
		a[rs].key+=a[cur].add;
		a[ls].min+=a[cur].add;
		a[rs].min+=a[cur].add;
		a[cur].add=0;
	}
	if(a[cur].rev>0)
	{
		a[ls].rev^=1;
		a[rs].rev^=1;
		a[cur].le=rs;
		a[cur].ri=ls;
		a[cur].rev=0;
	}
}

void update(int cur)
{
	int ls=a[cur].le,rs=a[cur].ri;
	a[cur].size=a[ls].size+a[rs].size+1;
	a[cur].min=a[cur].key;
	if(ls&&a[ls].min<a[cur].min)a[cur].min=a[ls].min;
	if(rs&&a[rs].min<a[cur].min)a[cur].min=a[rs].min;
}	

void leftrotate(int x)
{
	int y=a[x].ri,p=a[x].pre;
	a[x].ri=a[y].le;
	if(a[x].ri)a[a[x].ri].pre=x;
	a[y].le=x;
	a[x].pre=y;
	a[y].pre=p;
	if(!p)T=y;
	else
		a[p].ri==x?a[p].ri=y:a[p].le=y;
	update(x);
}

void rightrotate(int x)
{
	int y=a[x].le,p=a[x].pre;
	a[x].le=a[y].ri;
	if(a[x].le)a[a[x].le].pre=x;
	a[y].ri=x;
	a[x].pre=y;
	a[y].pre=p;
	if(!p)T=y;
	else
		a[p].ri==x?a[p].ri=y:a[p].le=y;
	update(x);
}

void splay(int x,int goal)
{
	int y,z;
	while(1)
	{
		if((y=a[x].pre)==goal)break;
		if((z=a[y].pre)==goal)
			a[y].ri==x?leftrotate(y):rightrotate(y);
		else
		{
			if(a[z].ri==y)
			{
				if(a[y].ri==x)
					leftrotate(z),leftrotate(y);
				else
					rightrotate(y),leftrotate(z);
			}
			else
			{
				if(a[y].le==x)
					rightrotate(z),rightrotate(y);
				else
					leftrotate(y),rightrotate(z);
			}
		}
	}
	update(x);
}

void rotateto(int k,int goal)
{
	int i=T;
	while(1)
	{
		pushdown(i);
		if(a[a[i].le].size+1==k)break;
		if(k<=a[a[i].le].size)i=a[i].le;
		else k-=a[a[i].le].size+1,i=a[i].ri;
	}
	splay(i,goal);
}

void newnode(int &cur,int v)
{
	cur=++node;
	a[cur].min=a[cur].key=v;
	a[cur].size=1;
	a[cur].le=a[cur].ri=a[cur].rev=a[cur].add=0;
}

void build(int &cur,int x,int y,int p)
{
	int mid=(x+y)>>1;
	newnode(cur,s[mid]);
	a[cur].pre=p;
	if(x==y)return;
	if(x<mid)build(a[cur].le,x,mid-1,cur);
	if(y>mid)build(a[cur].ri,mid+1,y,cur);
	update(cur);
}

void init(int n)
{
	int i;
	memset(s,0,sizeof(s));
	memset(a,0,sizeof(a));
	for(i=1;i<=n;i++)scanf("%d",&s[i]);
	T=node=0;
	build(T,0,n+1,0);
}

void Add(int x,int y,int z)
{
	int k;
	rotateto(x,0);rotateto(y+2,T);
	k=a[a[T].ri].le;
	a[k].add+=z;a[k].key+=z;a[k].min+=z;
}

void Reverse(int x,int y)
{
	int k;
	rotateto(x,0);rotateto(y+2,T);
	k=a[a[T].ri].le;
	a[k].rev^=1;
}

void Revolve(int x,int y,int z)
{
	int k=z%(y-x+1),t;
	if(k)
	{
		rotateto(x,0);rotateto(y-k+2,T);
		t=a[a[T].ri].le;
		a[a[T].ri].le=0;
		update(a[T].ri);update(T);
		rotateto(x+k,0);rotateto(x+k+1,T);
		a[a[T].ri].le=t;a[t].pre=a[T].ri;
		update(a[T].ri);update(T);
	}
}

void Insert(int x,int y)
{
	rotateto(x+1,0);rotateto(x+2,T);
	newnode(a[a[T].ri].le,y);
	a[a[a[T].ri].le].pre=a[T].ri;
	update(a[T].ri);update(T);
}

void Delete(int x)
{
	rotateto(x,0);rotateto(x+2,T);
	a[a[T].ri].le=0;
	update(a[T].ri);update(T);
}

void Min(int x,int y)
{
	rotateto(x,0);rotateto(y+2,T);
	printf("%d\n",a[a[a[T].ri].le].min);
}
```

##树链剖分

###树链剖分(边)

```cpp
vector<pair<int, int>> v[200001];//边及该边的编号
int w[200001];//边权
int n, cnt;
int father[200001], depth[200001], top[200001], id[200001];
int f[200001];//边在树状数组（线段树）中的位置
int tmp[200001];
int dfs1(int i, int fa)
{
	father[i] = fa;
	depth[i] = depth[fa] + 1;
	tmp[i] = -1;
	int ret = 0, maxSize = 0;
	for (unsigned int j = 0; j < v[i].size(); j++){
		int t = v[i][j].first;
		if (t == fa)continue;
		int size = dfs1(t, i);
		ret += size;
		if (size > maxSize){
			maxSize = size;
			tmp[i] = j;
		}
	}
	return ret + 1;
}
void dfs2(int i, int tp, int index)
{
	top[i] = tp;
	id[i] = cnt;
	f[index] = cnt++;
	if (tmp[i] != -1)
		dfs2(v[i][tmp[i]].first, tp, v[i][tmp[i]].second);
	for (unsigned int j = 0; j < v[i].size(); j++){
		int t = v[i][j].first;
		if (t != father[i] && j != tmp[i])
			dfs2(t, t, v[i][j].second);
	}
}
int queryTree(int s, int t)
{
	int ret = 0;
	int top1 = top[s], top2 = top[t];
	while (top1 != top2){
		if (depth[top1] < depth[top2]){
			ret += sum(id[t]) - sum(id[top2] - 1);
			t = father[top2]; top2 = top[t];
		}
		else{
			ret += sum(id[s]) - sum(id[top1] - 1);
			s = father[top1]; top1 = top[s];
		}
	}
	if (s != t){
		if (depth[s] > depth[t])swap(s, t);
		ret += sum(id[t]) - sum(id[s]);
	}
	return ret;
}
void init()
{
	cnt = 0;
	dfs1(1, 1);
	dfs2(1, 1, 0);
	for (int i = 1; i < n; i++)
		tree[f[i]] = w[i];
	build();
}
int main()
{
	int q, cur;
	scanf("%d%d%d", &n, &q, &cur);
	for (int i = 1; i < n; i++){
		int s, t, value;
		scanf("%d%d%d", &s, &t, &value);
		v[s].push_back(make_pair(t, i));
		v[t].push_back(make_pair(s, i));
		w[i] = value;
	}
	init();
	for (int i = 0; i < q; i++){
		int t, u, value;
		scanf("%d%d", &t, &u);
		if (t == 0){
			printf("%d\n", queryTree(cur, u));
			cur = u;
		}
		else{
			scanf("%d", &value);
			add(f[u], value - w[u]);
			w[u] = value;
		}
	}
	return 0;
}
```

###树链剖分(点)

```cpp
vector<int> v[100001];
int n, cnt, color;
int father[100001], depth[100001], top[100001], id[100001], son[100001];
struct Tree{
	int maxValue, maxId, delta;
	bool set;
}tree[1 << 18];
int treeLen;
int dfs1(int i, int fa)
{
	father[i] = fa;
	depth[i] = depth[fa] + 1;
	son[i] = 0;
	int ret = 0, maxSize = 0;
	for (unsigned int j = 0; j < v[i].size(); j++){
		int t = v[i][j];
		if (t == fa)continue;
		int size = dfs1(t, i);
		ret += size;
		if (size > maxSize){
			maxSize = size;
			son[i] = t;
		}
	}
	return ret + 1;
}
void dfs2(int i, int tp)
{
	top[i] = tp;
	id[i] = cnt++;
	if (son[i])dfs2(son[i], tp);
	for (unsigned int j = 0; j < v[i].size(); j++){
		int t = v[i][j];
		if (t != father[i] && t != son[i])dfs2(t, t);
	}
}
void init()
{
	cnt = 0; depth[1] = 0;
	dfs1(1, 1);
	dfs2(1, 1);
	for (treeLen = 1; treeLen < n; treeLen *= 2);
	memset(tree, 0, sizeof(Tree) * 2 * treeLen);
}
void pushDown(int i)
{
	if (tree[i].set){
		tree[2 * i].set = tree[2 * i + 1].set = true;
		tree[2 * i].delta = tree[2 * i + 1].delta = tree[i].delta;
		tree[i].delta = 0; tree[i].set = false;
	}
	else if (tree[i].delta){
		tree[2 * i].delta += tree[i].delta;
		tree[2 * i + 1].delta += tree[i].delta;
		tree[i].delta = 0;
	}
}
int queryL, queryR;
void addInternal(int i, int l, int len)
{
	if (queryL <= l && queryR >= l + len){
		tree[i].delta++;
		return;
	}
	len >>= 1; pushDown(i);
	int mid = l + len;
	if (mid > queryL)addInternal(2 * i, l, len);
	if (mid < queryR)addInternal(2 * i + 1, mid, len);
}
inline void addValue(int l, int r){
	queryL = l; queryR = r;
	addInternal(1, 0, treeLen);
}
int addTree(int s, int t)
{
	int ret = 0;
	int top1 = top[s], top2 = top[t];
	while (top1 != top2){
		if (depth[top1] < depth[top2]){
			addValue(id[top2], id[t] + 1);
			t = father[top2]; top2 = top[t];
		}
		else{
			addValue(id[top1], id[s] + 1);
			s = father[top1]; top1 = top[s];
		}
	}
	if (depth[s] > depth[t])swap(s, t);
	addValue(id[s], id[t] + 1);
	return ret;
}
void process(int i)
{
	if (tree[i].set){
		if (tree[i].delta > tree[i].maxValue){
			tree[i].maxValue = tree[i].delta;
			tree[i].maxId = color;
		}
		return;
	}
	pushDown(i);
	process(2 * i);
	process(2 * i + 1);
}
inline void pushDownColor(int i, int j){
	if (tree[j].maxValue < tree[i].maxValue
		|| (tree[j].maxValue == tree[i].maxValue && tree[i].maxId < tree[j].maxId)){
		tree[j].maxValue = tree[i].maxValue;
		tree[j].maxId = tree[i].maxId;
	}
}
void getAns(int i)
{
	if (i < treeLen){
		pushDownColor(i, 2 * i);
		pushDownColor(i, 2 * i + 1);
		getAns(2 * i);
		getAns(2 * i + 1);
	}
}
vector<pair<int, int>> z[100001];
int main()
{
	int m;
	while (scanf("%d%d", &n, &m) == 2 && n){
		for (int i = 1; i <= n; i++)
			v[i].clear();
		for (int i = 1; i < n; i++){
			int s, t;
			scanf("%d%d", &s, &t);
			v[s].push_back(t);
			v[t].push_back(s);
		}
		for (int i = 0; i < m; i++){
			int s, t, w;
			scanf("%d%d%d", &s, &t, &w);
			z[w].push_back(make_pair(s, t));
		}
		init();
		for (color = 1; color <= 100000; color++){
			tree[1].set = true; tree[1].delta = 0;
			for (unsigned int j = 0; j < z[color].size(); j++)
				addTree(z[color][j].first, z[color][j].second);
			process(1);
			z[color].clear();
		}
		getAns(1);
		for (int i = 1; i <= n; i++)
			printf("%d\n", tree[treeLen + id[i]].maxId);
	}
	return 0;
}
```

###Mergeable_Heap

```cpp
struct node
{
	int v,dis;
	node *l,*r;
}mem[maxn],*head[maxn];
int cnt;
node* merge(node* a,node* b)
{
	if (a==mem) return b;
	if (b==mem) return a;
	if (a->v<b->v) swap(a,b);
	a->r=merge(a->r,b);
	if (a->r->dis>a->l->dis) swap(a->l,a->r);
	if (a->r==mem) a->dis=0;
	else a->dis=a->r->dis+1;
	return a;
}
void init()
{
	mem[0].dis=-1;
	mem[0].l=mem[0].r=mem;
	for (int i=1;i<=n;i++) 
	{
		mem[i].l=mem[i].r=mem;
		head[i]=mem+i;
	}
}
//BZOJ 2809
int m;
queue<int> q;
int f[maxn],c[maxn],l[maxn],ind[maxn],ths[maxn];
long long cost[maxn],ans;
int main()
{
	scanf("%d%d",&n,&m);
	init();
	for (int i=1;i<=n;i++) 
	{
		scanf("%d%d%d",&f[i],&c[i],&l[i]);
		ind[f[i]]++;
	}
	for (int i=1;i<=n;i++)
	{
		if (ind[i]==0) q.push(i);
		ths[i]=1;cost[i]=c[i];
		mem[i].v=c[i];
	}
	while(!q.empty())
	{
		int now=q.front();q.pop();
		while(cost[now]>m)
		{
			cost[now]-=head[now]->v;
			head[now]=merge(head[now]->l,head[now]->r);
			ths[now]--;
		}
		ans=max(ans,1ll*l[now]*ths[now]);
		if (f[now]!=0)
		{
			head[f[now]]=merge(head[f[now]],head[now]);
			ths[f[now]]+=ths[now];
			cost[f[now]]+=cost[now];
			if ((--ind[f[now]])==0) q.push(f[now]);
		}
	}
	printf("%lld",ans);
	return 0;
}
```

***



***

***





# String

***

***

###Manacher

```cpp
gets(s1);l=0;
while(s1[0]!='E')
{
	l++;
	n=strlen(s1);
	s2[0]='$';k=0;
	for(i=0;i<n;i++)
	{
		s2[++k]='#';
		s2[++k]=s1[i];
	}
	s2[++k]='#';s2[++k]='\0';
	memset(p,0,sizeof(p));
	mx=0;id=0;
	for(i=1;s2[i]!='\0';i++)
	{
		p[i]=mx>i?min(p[2*id-i],mx-i):1;
		while(s2[i+p[i]]==s2[i-p[i]])p[i]++;
		if(i+p[i]>mx)
		{
			mx=i+p[i];id=i;
		}
	}
	mx=0;
	for(i=1;s2[i]!='\0';i++)if(p[i]-1>mx)mx=p[i]-1;
	printf("Case %d: %d\n",l,mx);
	gets(s1);
}
```

###KMP

```cpp
void get_nxt()
{
	int n=strlen(target_string);
	nxt[0]=0;nxt[1]=0;
	for (int i=1;i<n;i++)
	{
		int j=nxt[i];
		while(j&&target_string[i]!=target_string[j]) j=nxt[j];
		nxt[i+1]=target_string[i]==target_string[j]?j+1:0;
	}
}
int kmp()
{
	int n=strlen(origin_string);
	int m=strlen(target_string);
	int j=0,cnt=0;
	for (int i=0;i<n;i++)
	{
		while(j&&origin_string[i]!=target_string[j]) j=nxt[j];
		if (origin_string[i]==target_string[j]) j++;
		if (j==m) {cnt++;j=nxt[j];}
	}
	return cnt;
}
```

###字母表压缩版

```cpp
#include<cstdio>
#include<cstring>
#include<queue>
using namespace std;
#define LETTER 26
struct Trie{
	int num, next, fail;
}trie[1000000];
int cnt;
int pool[LETTER * 200000], poolEnd;
void init()
{
	cnt = 0;
	trie[0].num = 0;
	trie[0].next = -1;
	memset(pool, 0, sizeof(pool));
	poolEnd = 0;
}
inline int convert(char ch){ return ch - 'a'; }
inline bool oneBranch(int value){ return value < LETTER; }
inline int child(int i, int ch){
	if (oneBranch(trie[i].next))return trie[i].next == ch ? i + 1 : 0;
	return pool[trie[i].next + ch];
}
void insert(char *s)
{
	int pos = 0, i;
	for (i = 0; s[i]; i++){
		int t = trie[pos].next;
		if (oneBranch(t)){
			if (t == convert(s[i]))pos++;
			else{
				trie[pos].next = (poolEnd += LETTER);
				if (t != -1)pool[trie[pos].next + t] = pos + 1;
				break;
			}
		}
		else if (pool[t + convert(s[i])])
			pos = pool[t + convert(s[i])];
		else break;
	}
	if (s[i]){
		pool[trie[pos].next + convert(s[i])] = ++cnt;
		for (i++; s[i]; i++, cnt++){
			trie[cnt].num = 0;
			trie[cnt].next = convert(s[i]);
		}
		trie[cnt].num = 1;
		trie[cnt].next = -1;
	}
	else trie[pos].num++;
}
int getFailPoint(int father, int ch)
{
	while (father){
		father = trie[father].fail;
		int pos = child(father, ch);
		if (pos)return pos;
	}
	return 0;
}
void makeFail()
{
	queue<int> q; q.push(0);
	trie[0].fail = 0;
	while (!q.empty()){
		int t = q.front(); q.pop();
		if (oneBranch(trie[t].next)){
			if (trie[t].next != -1){
				trie[t + 1].fail = getFailPoint(t, trie[t].next);
				q.push(t + 1);
			}
		}
		else for (int i = 0; i < LETTER; i++){
			int cur = pool[trie[t].next + i];
			if (cur){
				trie[cur].fail=getFailPoint(t, i);
				q.push(cur);
			}
		}
	}
}
//统计匹配总次数，包括母串多次匹配同一模式串或多个模式串相同
int search(char *s)
{
	int ret = 0, cur = 0;
	for (int i = 0; s[i]; i++){
		int ch = convert(s[i]);
		for (; cur && !child(cur, ch); cur = trie[cur].fail);
		cur = child(cur, ch);
		for (int temp = cur; temp; temp = trie[temp].fail)
			ret += trie[temp].num;
	}
	return ret;
}
```

###回文树

```cpp
struct node
{
	int len,sum;
	node* fail,*next[26];
}mem[100005],*headf,*heads,*last;
int tot,now;
char s[100005];

void init()
{
	memset(mem,0,sizeof(mem));
	headf=mem;last=heads=mem+1;
	headf->fail=heads;heads->len=-1;
	tot=1;now=0;
}
void add(int x,int p)
{
	node* cur=last;
	for (;s[p-cur->len-1]!=s[p];cur=cur->fail);
	if (!cur->next[x])
	{
		node* ths=&mem[++tot];
		last=cur->next[x]=ths;
		ths->len=cur->len+2;
		if (cur==heads) ths->fail=headf;
		else
		{
			for (cur=cur->fail;s[p-cur->len-1]!=s[p];cur=cur->fail);
			ths->fail=cur->next[x];
		}
		ths->sum=ths->fail->sum+1;
	}
	else last=cur->next[x];
}
//HDU 5157
long long l[100005],r[100005];
int main()
{
	while(~scanf("%s",s))
	{
		int n=strlen(s);
		init();
		for (int i=0;i<n;i++) {add(s[i]-'a',i);l[i]=last->sum;}
		reverse(s,s+n);
		init();
		for (int i=0;i<n;i++) {add(s[i]-'a',i);r[i]=last->sum+r[i-1];}
		long long ans=0;
		for (int i=0;i<n-1;i++) ans+=l[i]*r[n-i-2];
		printf("%I64d\n",ans);
	}	
	return 0;
}
```



***

***





# Math

```cpp
void swap(double& p,double& q)
{
	double t;
	t=p;p=q;q=t;
}

struct Matrix
{
	double a[maxn][maxn];
	//1-n行表示第1-n个方程
	//每行第1-n个元素表示系数，第n+1个元素表示等号右边的常数
}q;

int ii,jj,nn;

LL det(LL a[][maxn], int n) {//求行列式值(整数版)
    int i, j, k, r;
    LL res = 1;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            while (a[j][i]) {
                LL f = a[i][i] / a[j][i];
                for (int k = i; k < n; k++) a[i][k] -= f * a[j][k];
                for (int k = i; k < n; k++) swap(a[i][k], a[j][k]);
                res = -res;
            }
        }
        if (a[i][i] == 0) return 0;
        res *= a[i][i];
    }
    return res < 0 ? -res : res;
}

double FF(double x)//需积分的函数，自行修改
{
	return 1.0;
}

double simpson(double x,double y)
{
	double z=x+(y-x)/2.0;
	return (y-x)/6.0*(FF(x)+FF(y)+4*FF(z));
}

double asr(double x,double y,double eeps,double A)//eeps为精度
{
	double z=x+(y-x)/2.0;
	double L=simpson(x,z);
	double R=simpson(z,y);
	if(fabs(L+R-A)<=15*eeps)return (L+R)+(L+R-A)/15.0;
	else return asr(x,z,eeps/2.0,L)+asr(z,y,eeps/2.0,R);
}

double simpson_zsx(double x,double y,double eeps)//自适应辛普森主函数
{
	return asr(x,y,eeps,simpson(x,y));
}

void gauss_eli(struct Matrix& p,int n)//高斯消元
{
	int i,j,k,r;
	for(i=1;i<=n;i++)
	{
		r=i;
		for(j=i+1;j<=n;j++)
			if(fabs(p.a[j][i])>fabs(p.a[r][i]))r=j;
		if(r!=i)for(j=1;j<=n+1;j++)swap(p.a[r][j],p.a[i][j]);
		for(k=1;k<=i-1;k++)
		{
			if(p.a[i][k]==0)continue;
			for(j=n+1;j>=k;j--)
				p.a[i][j]-=p.a[k][j]/p.a[k][k]*p.a[i][k];
		}
	}
	for(i=n;i>=1;i--)
	{
		for(j=i+1;j<=n;j++)
			p.a[i][n+1]-=p.a[j][n+1]*p.a[i][j];
		p.a[i][n+1]/=p.a[i][i];
	}
}

LL gcd(LL a,LL b)
{
	return b==0?a:gcd(b,a%b);
}

void tgcd(LL a,LL b,LL& d,LL& x,LL& y)//拓展欧几里德
{
	if(!b){d=a;x=1;y=0;}
	else{tgcd(b,a%b,d,y,x);y-=x*(a/b);}
}

LL pow_mod(LL a,LL p,LL n)//同余快速幂
{
	if(p==0)return 1;
	LL ans=pow_mod(a,p/2,n);
	ans=(ans*ans)%n;
	if(p%2==1)ans=(ans*a)%n;
	return ans;
}

int euler_phi(int n)//求欧拉函数
{
	int m=(int)sqrt(n+0.5);
	int ans=n;
	for(int i=2;i<=m;i++)
		if(n%i==0)
		{
			ans=ans/i*(i-1);
			while(n%i==0)n=n/i;
		}
	if(n>1)ans=ans/n*(n-1);
	return ans;
}

void phi_table(int n)//欧拉函数表
{
	memset(phi,0,n+1);
	phi[1]=1;
	for(int i=2;i<=n;i++)
	{
		if(phi[i])continue;
		for(int j=i;j<=n;j+=i)
		{
			if(!phi[j])phi[j]=j;
			phi[j]=phi[j]/i*(i-1);
		}
	}
}

LL inv(LL a,LL n)//a关于n的逆元
{
	LL d,x,y;
	tgcd(a,n,d,x,y);
	return d==1?(x+n)%n:-1;
}

LL china(int n,int* a,int* m)//中国剩余定理
{
	LL M=1,d,y,x=0;
	for(int i=0;i<n;i++)M*=m[i];
	for(int i=0;i<n;i++)
	{
		LL w=M/m[i];
		tgcd(m[i],w,d,d,y);
		x=(x+y*w*a[i])%M;
	}
	return (x+M)%M;
}

int log_mod(int a,int b,int n)//求解模方程a^x=b(mod n),n为素数,无解返回-1
{
	int m,v,e=1,i;
	m=(int)sqrt(n+0.5);
	v=inv(pow_mod(a,m,n),n);
	map<int,int> x;
	x[1]=0;
	for(i=1;i<m;i++)
	{
		e=(e*a)%n;
		if(!x.count(e))x[e]=i;
	}
	for(i=0;i<m;i++)
	{
		if(x.count(b))return (i*m+x[b]);
		b=(b*v)%n;
	}
	return -1;
}
```

###矩阵快速幂

```cpp
struct mat
{
    int n;
    LL num[105][105];

    void init0(int t)
    {
        n=t;
        for(int i=0;i<=n;i++)
            for(int j=0;j<=n;j++)
                num[i][j]=0;
    }

    void init1(int t)
    {
        n=t;
        for(int i=0;i<=n;i++)
            for(int j=0;j<=n;j++)
                if(i!=j)num[i][j]=0;else num[i][j]=1;
    }

    mat operator * (const struct mat p)const
    {
        struct mat ans;
        ans.init0(n);
        for(int i=1;i<=n;i++)
            for(int j=1;j<=n;j++)
                for(int k=1;k<=n;k++)
                    ans.num[i][j]=(ans.num[i][j]+num[i][k]*p.num[k][j])%mod;
        //printf("??");ans.testprint();
        return ans;
    }       

    mat operator ^ (int t)const
    {
        struct mat ans,now;
        ans.init1(n);
        now.n=n;
        for(int i=0;i<=n;i++)
            for(int j=0;j<=n;j++)
                now.num[i][j]=num[i][j];
        while(t>0)
        {
            if(t&1)ans=ans*now;
            now=now*now;
            t>>=1;
        }
        return ans;
    }

};
```

***

## 傅里叶变换

###FFT

```cpp
complex<double> epsilon[maxn];
complex<double> arti_epsilon[maxn];
complex<double> a[maxn],b[maxn],c[maxn],temp[maxn];

int n1,n2,m;

void init_epsilon(int n)
{
	for(int i=0;i!=n;i++)
	{
		epsilon[i]=complex<double>(cos(2.0*pi*i/n),sin(2.0*pi*i/n));
		arti_epsilon[i]=conj(epsilon[i]);
	}
}

int calc(int t)
{
	int j=0;
	while((1<<j)<=t)j++;
	return 1<<j;
}

void DFT(int n,complex<double>*  buffer,int offset,int step,complex<double>* epsilon)
{
	if(n==1)return;
	int m=n>>1;
	DFT(m,buffer,offset,step<<1,epsilon);
	DFT(m,buffer,offset+step,step<<1,epsilon);
	for(int k=0;k!=m;k++)
	{
		int pos=2*step*k;
		temp[k]=buffer[pos+offset]+epsilon[k*step]*buffer[pos+offset+step];
		temp[k+m]=buffer[pos+offset]-epsilon[k*step]*buffer[pos+offset+step];
	}
	for(int i=0;i!=n;i++)buffer[i*step+offset]=temp[i];
}

//IDFT 将 DFT 中的epsilon改为arti_epsilon即可
	
void FFT(int m,complex<double>* a,complex<double>* b,complex<double>* c)
{
	init_epsilon(m);
	DFT(m,a,0,1,epsilon);
	DFT(m,b,0,1,epsilon);
	for(int i=0;i<=m;i++)c[i]=a[i]*b[i];
	IDFT(m,c,0,1,epsilon);
	double mm=m;
	for(int i=0;i<=m;i++)c[i]/=mm;
}

int init()//n1,n2表示多项式次数
{
	double x,y;
	scanf("%d%d",&n1,&n2);
	memset(a,0,sizeof(a));
	memset(b,0,sizeof(b));
	for(int i=0;i<=n1;i++)
	{
		scanf("%lf %lf",&x,&y);
		a[i].real(x);
		a[i].imag(y);
	}
	for(int i=0;i<=n2;i++)
	{
		scanf("%lf %lf",&x,&y);
		b[i].real(x);
		b[i].imag(y);
	}
	m=calc(n1+n2);
	return m;
}

void print()
{
	for(int i=0;i<m;i++)printf("%lf %lf\n",real(c[i]),imag(c[i]));
}

int main()
{
	m=init();
	FFT(m,a,b,c);
	print();
}
```

###NTT&CRT

```cpp
int len, bit;
int MOD, w[2][32];
inline int add(int a, int b){
	return a + b - (a + b >= MOD ? MOD : 0);
}
inline int sub(int a, int b){
	return a - b + (a - b < 0 ? MOD : 0);
}
inline int mul(int a, int b){
	return (long long)a * b % MOD;
}
int power(int a, int b){
	int ret = 1;
	for (int t = a; b; b >>= 1){
		if (b & 1)ret = mul(ret, t);
		t = mul(t, t);
	}
	return ret;
}
int cal_root(int mod)
{
	for (int i = 2;; i++){
		if (power(i, (mod - 1) / 2) == mod - 1)
			return i;
	}
}
void fft_init(int n, int mod)
{
	MOD = mod;
	bit = (int)log2(n - 0.5) + 2;
	len = 1 << bit;
	w[0][0] = power(cal_root(mod), (mod - 1) / len);
	int i;
	for (i = 1; i < bit; i++)
		w[0][i] = mul(w[0][i - 1], w[0][i - 1]);
	i--;
	w[1][i] = w[0][i];
	for (i--; i >= 0; i--)
		w[1][i] = mul(w[1][i + 1], w[0][i]);
}
void bitReverse(int a[]) {
	for (int i = 1, j = len / 2; i < len - 1; i++) {
		if (i < j) swap(a[i], a[j]);
		int k = len / 2;
		while (j >= k) { j -= k; k >>= 1; }
		if (j < k) j += k;
	}
}
void fft_main(int a[], bool reverse)
{
	bitReverse(a);
	for (int i = 1, s = 1; s < len; i++, s <<= 1){
		int step = w[reverse][bit - i];
		for (int j = 0; j < len; j += 2 * s){
			int cur = 1;
			for (int k = j; k < j + s; k++){
				int u = a[k], t = mul(cur, a[k + s]);
				a[k] = add(u, t);
				a[k + s] = sub(u, t);
				cur = mul(cur, step);
			}
		}
	}
	if (reverse){
		int t = power(len, MOD - 2);
		for (int i = 0; i < len; i++)
			a[i] = mul(a[i], t);
	}
}
//确保数组中的数小于mod(mod<2^30)，数组需留足2^(logn向上取整+1)的空间，后面填充0
//并且mod为形如m*2^k+1的素数，2^k>=2*n
void fft(int a[], int b[], int n, int mod)
{
	fft_init(n, mod);
	fft_main(a, 0); fft_main(b, 0);
	for (int i = 0; i < len; i++)
		a[i] = mul(a[i], b[i]);
	fft_main(a, 1);
}
//确保mod两两互质，retmod任意
void chineseRemainder(const int mod[], int *a[], int ret[], int num, int n, int retMod)
{
	int kk[30], mulMod[30][30], mulModr[30], mulretMod[30];
	for (int i = 0; i < num; i++){
		MOD = mod[i]; mulMod[i][0] = 1;
		for (int j = 1; j <= i; j++)
			mulMod[i][j] = mul(mulMod[i][j - 1], mod[j - 1]);
		mulModr[i] = power(mulMod[i][i], MOD - 2);
	}
	mulretMod[0] = 1; MOD = retMod;
	for (int i = 1; i < num; i++)
		mulretMod[i] = mul(mulretMod[i - 1], mod[i - 1]);
	for (int i = 0; i < n; i++){
		for (int j = 1; j < num; j++){
			MOD = mod[j];
			int sum = a[0][i] % MOD;
			for (int k = 1; k < j; k++)
				sum = add(sum, mul(mulMod[j][k], kk[k]));
			kk[j] = mul(sub(a[j][i] % MOD, sum), mulModr[j]);
		}
		MOD = retMod;
		ret[i] = a[0][i] % MOD;
		for (int j = 1; j < num; j++)
			ret[i] = add(ret[i], mul(kk[j] % MOD, mulretMod[j]));
	}
}
//附满足条件大整数：167772161, 469762049, 754974721
```

***

## 决战素数！

###Miller-Rabin&&Pollard

```cpp
const UInt base1[] = { 2, 7, 61, 0 };
const UInt base2[] = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022, 0 };
const UInt prime[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53 };
template <typename T>
inline T add(T a, T b, T mod){
	return a + b - (a + b >= mod ? mod : 0);
}
inline UInt mul(UInt a, UInt b, UInt mod){
	return (ULL)a * b % mod;
}
ULL mul(ULL a, ULL b, ULL mod){
	ULL ret = 0;
	for (ULL t = a; b; b >>= 1){
		if (b & 1)ret = add(ret, t, mod);
		t <<= 1;
		if (t >= mod)t -= mod;
	}
	return ret;
}
template <typename T>
T power(T a, T b, T mod){
	T ret = 1;
	for (T t = a; b; b >>= 1){
		if (b & 1)ret = mul(ret, t, mod);
		t = mul(t, t, mod);
	}
	return ret;
}
//n为小于2^63的非1奇数，正确性100%
template <typename T>
bool millerRabin(T n)
{
	int s = 0;
	T r = n;
	for (r--; !(r & 1); r >>= 1)s++;
	for (const UInt *base = typeid(T) == typeid(UInt) ? base1 : base2; *base; base++){
		T t = power(*base % n, r, n);
		if (t == 0 || t == 1 || t == n - 1)continue;
		for (int j = 1; j < s; j++){
			t = mul(t, t, n);
			if (t == 1)return false;
			if (t == n - 1)break;
		}
		if (t != n - 1)return false;
	}
	return true;
}
template <typename T>
bool checkPrime(T n)
{
	if (n == 1)return false;
	for (int i = 0; i < sizeof(prime) / sizeof(int); i++){
		if (n % prime[i] == 0)return n == prime[i];
	}
	return millerRabin(n);
}
template <typename T>
T gcd(T x, T y){
	return y ? gcd(y, x % y) : x;
}
template <typename T>
T pollard(T n)
{
	if (millerRabin(n))return n;
	while (1){
		T x = rand() % n, y = x, c = rand() % (n - 1) + 1;
		for (UInt i = 1, j = 2;; i++){
			if (i == j){ j *= 2; y = x; }
			x = add(mul(x, x, n), c, n);
			T d = gcd(x - y + n, n);
			if (d != 1){
				if (d != n)return d;
				break;
			}
		}
	}
}
ULL factor[64];
int factorNum;
void calFactorInternal(ULL n)
{
	ULL d;
	d = n >> 32 ? pollard(n) : pollard((UInt)n);
	if (d == n){ factor[factorNum++] = d; return; }
	calFactorInternal(d);
	calFactorInternal(n / d);
}
void calFactor(ULL n)
{
	factorNum = 0;
	for (int i = 0; i < sizeof(prime) / sizeof(int); i++){
		while (n % prime[i] == 0){
			n /= prime[i];
			factor[factorNum++] = prime[i];
		}
	}
	if (n != 1)calFactorInternal(n);
	sort(factor, factor + factorNum);
}
```

###Euler筛法

```cpp
void calPrime()
{
	for (int i = 2; i < MAXN; i++){
		if (!minFactor[i]){
			prime[primeNum++] = i;
			minFactor[i] = primeNum;
		}
		for (int j = 1; j <= minFactor[i]; j++){
			int t = i * prime[j - 1];
			if (t >= MAXN)break;
			minFactor[t] = j;
		}
	}
}
void calPhi()
{
	phi[1] = 1;
	for (int i = 2; i < MAXN; i++){
		if (!minFactor[i]){
			prime[primeNum++] = i;
			minFactor[i] = primeNum;
			phi[i] = i - 1;
		}
		for (int j = 1;; j++){
			int t = i * prime[j - 1];
			if (t >= MAXN)break;
			minFactor[t] = j;
			if (j == minFactor[i]){
				phi[t] = phi[i] * prime[j - 1];
				break;
			}
			phi[t] = phi[i] * (prime[j - 1] - 1);
		}
	}
}
```

