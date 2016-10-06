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
