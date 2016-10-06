int n;						//µãÊý
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
	freopen("in.txt","r",stdin);
	freopen("out.txt","w",stdout);
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

