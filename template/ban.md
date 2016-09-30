# 字符串

### KMP

```cpp
	next[0]=-1;i=0;j=-1;
	while(i<nm)
	{
		if(j==-1||m[i]==m[j])
		{
			++i;++j;
			if(m[i]!=m[j])next[i]=j;else next[i]=next[j];
		}else j=next[j];                      
	};
	for(i=0,j=0;i<ns;i++)
	{
		while(s[i]!=m[j]&&j>0)j=next[j-1];
		if(s[i]==m[j])j++;
		if(j>=nm){printf("%d\n",i-j+1);j=next[j-1];}
	}
```

### ACMachine

```cpp
struct tree
{
	int fail,num,fat;
	int hi;
	int son[27];
}a[500005];
char s[55];
char ss[1000005];
int i,j,m,n,ans;
int _,__;
queue<int> q;

void mt(char *x)
{
	int ii,o,ll;
	o=1;
	ll=strlen(x);
	for(ii=0;ii<ll;ii++)
	{
		if(a[o].son[x[ii]-96]>0)o=a[o].son[x[ii]-96];
		else
		{
			m++;a[m].fat=o;
			a[o].son[x[ii]-96]=m;
			a[m].num=x[ii]-96;
			o=m;
		}
		if(ii>=ll-1)a[o].hi++;
	}
}

void ACM()
{
	while(!q.empty())
	{
		int r=q.front();
		for(int ii=1;ii<=26;ii++)
			if(a[r].son[ii]>0)q.push(a[r].son[ii]);
		if(r>1)
		{
			if(a[r].fat==1)a[r].fail=1;
			else
			{
				int t=a[r].fat;t=a[t].fail;
				while(t>1&&a[t].son[a[r].num]==0)t=a[t].fail;
				if(a[t].son[a[r].num]>0)a[r].fail=a[t].son[a[r].num];
				else a[r].fail=t;
			}
		}
		q.pop();
	}
	return;
}

void mat(char *s)
{
	int t,ii,ll,o;
	ii=0;ll=strlen(s);
	while(ii<ll)
	{
		while((ii<ll)&&(a[1].son[s[ii]-96]==0))ii++;
		if(ii>=ll)break;
		o=a[1].son[s[ii]-96];
		while(o>1)
		{
			if(a[o].hi>0)
			{
				ans+=a[o].hi;
				a[o].hi=0;
				t=a[o].fail;
				while(a[t].num==s[ii]-96)
				{
					if(a[t].hi>0){ans+=a[t].hi;a[t].hi=0;}
					t=a[t].fail;
				}
			}
			ii++;
			if(ii>=ll)break;
			if(a[o].son[s[ii]-96]>0)o=a[o].son[s[ii]-96];
			else
			{
				o=a[o].fail;
				while(o>1&&a[o].son[s[ii]-96]==0)o=a[o].fail;
				if(o>1)o=a[o].son[s[ii]-96];
			}
		}
	}
}
		m=1;memset(a,0,sizeof(a));
		if(!q.empty())while(!q.empty())q.pop();
		a[m].fail=a[m].num=a[m].fat=0;a[m].hi=0;
		scanf("%d\n",&n);
		for(i=1;i<=n;i++)
		{
			scanf("%s\n",s);
			mt(s);
		}
		q.push(1);ans=0;
		a[1].fail=0;ACM();
		scanf("%s\n",ss);
		mat(ss);
		printf("%d\n",ans);
```

### Manacher

```cpp
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
```

# 图论

### Dijkstra(nlogn)

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

### Dinic

```cpp
struct edge
{
	int n,v;
	edge* next;
}mem[2000005],*head[2000],*store[505];
queue<int> q;
int cnt,n;
int dis[2000];
char ch;
void add_edge(int s,int t,int v)
{
	mem[cnt].n=t;mem[cnt].v=v;
	mem[cnt].next=head[s];head[s]=mem+cnt;
	cnt++;
	mem[cnt].n=s;mem[cnt].v=0;
	mem[cnt].next=head[t];head[t]=mem+cnt;
	cnt++;
}
bool bfs(int x)
{
	for (int i=0;i<=n;i++) dis[i]=INF;
	q.push(x);dis[x]=0;
	while(!q.empty())
	{
		for (edge *it=head[q.front()];it;it=it->next)
		if (it->v&&dis[it->n]==INF)
		{
			dis[it->n]=dis[q.front()]+1;
			q.push(it->n);
		}
		q.pop();
	}
	return (dis[n]!=INF);
}
int dinic(int x,int f)
{
	if (x==n) return f;
	for (edge* it=head[x];it;it=it->next)
	if (dis[it->n]==dis[x]+1&&it->v)
	{
		int nowf=dinic(it->n,min(f,it->v));
		if (nowf)
		{
			it->v-=nowf;
			mem[(it-mem)^1].v+=nowf;
			return nowf;
		}
	}
	return 0;
}
int DINIC(int x)
{
	int ths=0;
	while(bfs(x))
	{
		while(1)
		{
			int tmp=dinic(x,INF);
			if (!tmp) break;
			ths+=t;
		}
	}
	return ths;
}
//flow=DINIC(source);
```

# 数学

### 数论

```cpp
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

### FFT

```cpp
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

void IDFT(int n,complex<double>*  buffer,int offset,int step,complex<double>* epsilon)
{
	if(n==1)return;
	int m=n>>1;
	IDFT(m,buffer,offset,step<<1,epsilon);
	IDFT(m,buffer,offset+step,step<<1,epsilon);
	for(int k=0;k!=m;k++)
	{
		int pos=2*step*k;
		temp[k]=buffer[pos+offset]+arti_epsilon[k*step]*buffer[pos+offset+step];
		temp[k+m]=buffer[pos+offset]-arti_epsilon[k*step]*buffer[pos+offset+step];
	}
	for(int i=0;i!=n;i++)buffer[i*step+offset]=temp[i];
}
	
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
```

# 数据结构

### 区间修改线段树

```c++
struct tree
{
	int mi,ls,rs,ll,rr,add;
}a[3*maxn];

int c[maxn];
int i,j,k,l,m,n,t,T;

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



### Splay(Superemo)

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

### 笛卡尔树

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

### 堆

```cpp
void build()
{
	if ((n & 1) == 0 && a[n >> 1] > a[n])
		swap(a[n >> 1], a[n]);
	for (int i = (n - 1) / 2; i > 0; i--){
		if (a[2 * i] < a[i])swap(a[i], a[2 * i]);
		if (a[2 * i + 1] < a[i])swap(a[i], a[2 * i + 1]);
	}
}
void push(int value)
{
	a[++n] = value;
	for (int i = n, j = i >> 1; j; i = j, j >>= 1){
		if (a[i] < a[j])swap(a[i], a[j]);
	}
}
void pop()
{
	a[1] = a[n--];
	for (int i = 1, j = 2; j <= n; i = j, j *= 2){
		if (j != n && a[2 * i + 1] < a[2 * i])j++;
		if (a[j] < a[i])swap(a[i], a[j]);
	}
}
```

### 树状数组

```cpp
void build(int tree[], int m)
{
	for (int i = 2; i <= m; i *= 2){
		for (int j = i; j <= m; j += i)
			tree[j] += tree[j - i / 2];
	}
}
int sum(int tree[], int pos)
{
	int ret = 0;
	while (pos){
		ret += tree[pos];
		pos -= pos&-pos;
	}
	return ret;
}
void add(int tree[], int pos, int value, int m)
{
	while (pos <= m){
		tree[pos] += value;
		pos += pos&-pos;
	}
}
```

# 计算几何

```cpp
const double eps=1e-10;
const double pi=3.1415926535897932384626433832795;
const double eln=2.718281828459045235360287471352;
#include<bits/stdc++.h>
using namespace std;
#define For(i,n) for(int i=1;i<=n;i++)
#define Fork(i,k,n) for(int i=k;i<=n;i++)
#define Rep(i,n) for(int i=0;i<n;i++)
#define ForD(i,n) for(int i=n;i;i--)
#define ForkD(i,k,n) for(int i=n;i>=k;i--)
#define RepD(i,n) for(int i=n;i>=0;i--)
#define Forp(x) for(int p=Pre[x];p;p=Next[p])
#define Forpiter(x) for (int &p=iter[x];p;p=Next[p])
#define Lson (o<<1)
#define Rson ((o<<1)+1)
#define MEM(a) memset(a,0,sizeof(a));
#define MEMI(a) memset(a,127,sizeof(a));
#define MEMi(a) memset(a,128,sizeof(a));
#define INF (2139062143)
#define F (100000007)
#define pb push_back
#define mp make_pair
#define fi first
#define se second
#define vi vector<int>
#define pi pair<int,int>
#define SI(a) ((a).size())
typedef long long ll;
typedef long double ld;
typedef unsigned long long ull;
ll mul(ll a,ll b) {return (a*b)%F;}
ll add(ll a,ll b) {return (a+b)%F;}
ll sub(ll a,ll b) {return (a-b+llabs(a-b)/F*F+F)%F;}
void upd(ll &a,ll b){a=(a%F+b%F)%F;}
int read()
{
	int x=0,f=1; char ch=getchar();
	while(!isdigit(ch)) {if (ch=='-') f=-1;ch=getchar();}
	while(isdigit(ch)) {x=x*10+ch-'0';ch=getchar();}
	return x*f;
}
ll sqr(ll a){return a*a;}
ld sqr(ld a){return a*a;}
double sqr(double a){return a*a;}

ld PI=3.141592653589793238462643383;
class P{
public:
		double x,y;
		P(double x=0,double y=0):x(x),y(y){}
		friend long double dis2(P A,P B){return sqr(A.x-B.x)+sqr(A.y-B.y);}
		friend long double Dot(P A,P B){return A.x*B.x+A.y*B.y;}
		friend long double Length(P A){return sqrt(Dot(A,A));}
		friend long double Angle(P A,P B){return acos(Dot(A,B)/Length(A)/Length(B));}

		friend P operator- (P A,P B) {return P(A.x-B.x,A.y-B.y);}
		friend P operator+ (P A,P B) {return P(A.x+B.x,A.y+B.y);}
		friend P operator* (P A,double p) {return P(A.x*p,A.y*p);}
		friend P operator/ (P A,double p) {return P(A.x/p,A.y/p);}
		friend bool operator< (const P& a,const P& b) {return a.x<b.x||(a.x==b.x&&a.y<b.y);}

};
const double eps=1e-10;
int dcmp(double x){
	if (fabs(x)<eps) return 0;else return x<0?-1:1;
}
bool operator==(const P& a,const P& b){
	return dcmp(a.x-b.x)==0&&dcmp(a.y-b.y)==0;
}
typedef P V;
double Cross(V A,V B) {return A.x*B.y-A.y*B.x;}
double Area2(P A,P B,P C) {return Cross(B-A,C-A);}
V Rotate(V A,double rad){
	return V(A.x*cos(rad)-A.y*sin(rad),A.x*sin(rad)+A.y*cos(rad));
}
V Normal(V A){
	double L=Length(A);
	return V(-A.y/L,A.x/L);
}

P GetLineIntersection(P p,V v,P Q,V w){
	V u=p-Q;
	double t=Cross(w,u)/Cross(v,w);
	return p+v*t;
}
P GetLineIntersectionB(P p,V v,P Q,V w){
	return GetLineIntersection(p,v-p,Q,w-Q);
}

double DistanceToLine(P p,P A,P B){
	V v1=B-A,v2=p-A;
	return fabs(Cross(v1,v2))/Length(v1);
}
double DistanceToSegment(P p,P A,P B){
	if (A==B) return Length(p-A);
	V v1=B-A,v2=p-A,v3=p-B;
	if (dcmp(Dot(v1,v2))<0) return Length(v2);
	else if (dcmp(Dot(v1,v3))>0) return Length(v3);
	else return fabs(Cross(v1,v2))/Length(v1);
}
P GetLineProjection(P p,P A,P B){
	V v=B-A;
	return A+v*(Dot(v,p-A)/Dot(v,v));
}
bool SegmentProperIntersection(P a1,P a2,P b1,P b2){
	double c1=Cross(a2-a1,b1-a1),c2=Cross(a2-a1,b2-a1),c3=Cross(b2-b1,a1-b1),c4=Cross(b2-b1,a2-b1);
	return dcmp(c1)*dcmp(c2)<0&&dcmp(c3)*dcmp(c4)<0;
}
bool OnSegment(P p,P a1,P a2){
	return (dcmp(Cross(a1-p,a2-p))==0)&&(dcmp(Dot(a1-p,a2-p))<0);
}
double PolygonArea(P *p,int n){
	double area=0;
	For(i,n-2) area+=Cross(p[i]-p[0],p[i+1]-p[0]);
		return area/2;
}
P read_point(){
	P a;
	scanf("%lf%lf",&a.x,&a.y);
	return a;
}
struct C{
	P c;
	double r,x,y;
	C(P c,double r):c(c),r(r),x(c.x),y(c.y){}
	P point (double a){
		return P(c.x+cos(a)*r,c.y+sin(a)*r);	
	}
};

struct Line{
	P p;
	V v;
	double ang;
	Line(){}
	Line(P p,V v):p(p),v(v) {ang=atan2(v.y,v.x);}
	bool operator<(const Line & L) const{
		return ang<L.ang;
	}
	P point(double a){
		return p+v*a;
	}
};
int getLineCircleIntersection(Line L,C cir,double &t1,double &t2,vector<P> & sol){
	if (dcmp(DistanceToLine(cir.c,L.p,L.p+L.v)-cir.r)==0){
		P A=GetLineProjection(cir.c,L.p,L.p+L.v);
		sol.pb(A);
		t1=(A-L.p).x/L.v.x;
		return 1;
	}
	double a=L.v.x,b=L.p.x-cir.c.x,c=L.v.y,d=L.p.y-cir.c.y;
	double e=a*a+c*c,f=2*(a*b+c*d),g=b*b+d*d-cir.r*cir.r;
	double delta=f*f-4*e*g;
	if (dcmp(delta)<0) return 0;
	else if (dcmp(delta)==0){
		t1=t2=-f/(2*e);sol.pb(L.point(t1));
		return 1;
	}
	t1=(-f-sqrt(delta))/(2*e);sol.pb(L.point(t1));
	t2=(-f+sqrt(delta))/(2*e);sol.pb(L.point(t2));
	return 2;
}
double angle(V v){return atan2(v.y,v.x);}
int getCircleCircleIntersection(C C1,C C2,vector<P>& sol){
	double d=Length(C1.c-C2.c);
	if (dcmp(d)==0){
		if (dcmp(C1.r-C2.r)==0) return -1;
		return 0;
	}
	if (dcmp(C1.r+C2.r-d)<0) return 0;
	if (dcmp(fabs(C1.r-C2.r)-d)>0) return 0;

	double a=angle(C2.c-C1.c);
	double da=acos((C1.r*C1.r+d*d-C2.r*C2.r)/(2*C1.r*d));
	P p1=C1.point(a-da),p2=C1.point(a+da);
	sol.pb(p1);
	if (p1==p2) return 1;
	sol.pb(p2);
	return 2;
}
int getTangents(P p,C c,V* v){
	V u=c.c-p;
	double dist=Length(u);
	if (dist<c.r) return 0;
	else if (dcmp(dist-c.r)==0){
		v[0]=Rotate(u,PI/2);
		return 1;
	}else{
		double ang=asin(c.r/dist);
		v[0]=Rotate(u,-ang);
		v[1]=Rotate(u,ang);
		return 2;
	}
}
int getTangents(C A,C B,P* a,P* b){
	int cnt=0;
	if (A.r<B.r) {swap(A,B),swap(a,b);}
	int d2=(A.c.x-B.c.x)*(A.c.x-B.c.x)+(A.c.y-B.c.y)*(A.c.y-B.c.y);
	int rdiff=A.r-B.r;
	int rsum=A.r+B.r;
	if (d2<rdiff*rdiff) return 0;
	double base=atan2(B.y-A.y,B.x-A.x);
	if (d2==0&&A.r==B.r) return -1;
	if (d2==rdiff*rdiff){
		a[cnt]=A.point(base);b[cnt]=B.point(base);++cnt;
		return 1;
	}
	double ang=acos((A.r-B.r)/sqrt(d2));
	a[cnt]=A.point(base+ang);b[cnt]=B.point(base+ang);++cnt;
	a[cnt]=A.point(base-ang);b[cnt]=B.point(base-ang);++cnt;
	if (d2==rsum*rsum){
		a[cnt]=A.point(base);b[cnt]=B.point(PI+base);++cnt;
	}
	else if(d2>rsum*rsum){
		double ang=acos((A.r+B.r)/sqrt(d2));
		a[cnt]=A.point(base+ang);b[cnt]=B.point(PI+base+ang);++cnt;
		a[cnt]=A.point(base-ang);b[cnt]=B.point(PI+base-ang);++cnt;
	}
	return cnt;
}
C CircumscribedCircle(P p1,P p2,P p3){
	double Bx=p2.x-p1.x,By=p2.y-p1.y;
	double Cx=p3.x-p1.x,Cy=p3.y-p1.y;
	double D=2*(Bx*Cy-By*Cx);
	double cx=(Cy*(Bx*Bx+By*By)-By*(Cx*Cx+Cy*Cy))/D+p1.x;
	double cy=(Bx*(Cx*Cx+Cy*Cy)-Cx*(Bx*Bx+By*By))/D+p1.x;
	P p=P(cx,cy);
	return C(p,Length(p1-p));
}
C InscribedCircle(P p1,P p2,P p3){
	double a=Length(p2-p3);
	double b=Length(p3-p1);
	double c=Length(p1-p2);
	P p=(p1*a+p2*b+p3*c)/(a+b+c);
	return C(p,DistanceToLine(p,p1,p2));
}
double torad(double deg){
	return deg/180*acos(-1);
}
double radToPositive(double rad)
{
	if (dcmp(rad)<0) rad=ceil(-rad/PI)*PI+rad;
	if (dcmp(rad-PI)>=0) rad-=floor(rad/PI)*PI;
	return rad;
}
double todeg(double rad){
	return rad*180/acos(-1);
}
void get_coord(double R,double lat,double lng,double &x,double &y,double &z){
	lat=torad(lat);
	lng=torad(lng);
	x=R*cos(lat)*cos(lng);
	y=R*cos(lat)*sin(lng);
	z=R*sin(lat);
}
void print(double a){
	printf("%.6lf",a);
}
void print(P p){
	printf("(%.6lf,%.6lf)",p.x,p.y);
}
template<class T>
void print(vector<T> v){
	sort(v.begin(),v.end());
	putchar('[');
	int n=v.size();
	Rep(i,n){
		print(v[i]);
		if (i<n-1) putchar(',');
	}
	puts("]");
}
Line LineTranslation(Line l,V v){
	l.p=l.p+v;
	return l;
}
bool Convexpolygon(P &A,P &B,P &C,P &D){
	if (SegmentProperIntersection(A,C,B,D)) return 1;
	swap(B,C);
	if (SegmentProperIntersection(A,C,B,D)) return 1;
	swap(D,C);
	if (SegmentProperIntersection(A,C,B,D)) return 1;
	return 0;
}
bool IsParallel(P A,P B,P C,P D){
	return dcmp(Cross(B-A,D-C))==0;
}
bool IsPerpendicular(V A,V B){
	return dcmp(Dot(A,B))==0;
}
bool IsTrapezium(P A,P B,P C,P D){
	return IsParallel(A,B,C,D)^IsParallel(B,C,A,D);
}
bool IsParallelogram(P A,P B,P C,P D){
	return IsParallel(A,B,C,D)&&IsParallel(B,C,A,D);
}
bool IsRhombus(P A,P B,P C,P D){
	return IsParallelogram(A,B,C,D)&&dcmp(Length(B-A)-Length(C-B))==0;
}
bool IsRectangle(P A,P B,P C,P D){
	return IsParallelogram(A,B,C,D)&&IsPerpendicular(B-A,D-A);
}
bool IsSquare(P A,P B,P C,P D){
	return IsParallelogram(A,B,C,D)&&IsPerpendicular(B-A,D-A)&&dcmp(Length(B-A)-Length(C-B))==0;
}
double ArcDis(double chord,double r){
	return 2*asin(chord/2/r)*r;
}
typedef vector<P> Polygon;
int isPointInPolygon(P p,Polygon poly){
	int wn=0;
	int n=poly.size();
	Rep(i,n){
		if (OnSegment(p,poly[i],poly[(i+1)%n])) return -1;
		int k=dcmp(Cross(poly[(i+1)%n]-poly[i],p-poly[i]));
		int d1=dcmp(poly[i].y-p.y);
		int d2=dcmp(poly[(i+1)%n].y-p.y);
		if (k>0&&d1<=0&&d2>0) wn++;
		if (k<0&&d2<=0&&d1>0) wn--;
	}
	if (wn!=0) return 1;
	return 0;
}

int ConvexHull(P *p,int n,P *ch){
	sort(p,p+n);
	int m=0;
	Rep(i,n){
		while(m>1&&Cross(ch[m-1]-ch[m-2],p[i]-ch[m-2])<=0) m--;
		ch[m++]=p[i];
	}
	int k=m;
	RepD(i,n-2){
		while(m>k&&Cross(ch[m-1]-ch[m-2],p[i]-ch[m-2])<=0) m--;
		ch[m++]=p[i];
	}
	if (n>1) m--;
	return m;
}
void Two_pointFormToGeneralForm(P A,P B,double &a,double &b,double &c){
	a=A.y-B.y;
	b=B.x-A.x;
	c=Cross(A,B);
}
Polygon CutPolygon(Polygon poly,P A,P B){
	Polygon newpoly;
	int n=poly.size();
	Rep(i,n){
		P C=poly[i];
		P D=poly[(i+1)%n];
		if (dcmp(Cross(B-A,C-A))>=0) newpoly.pb(C);
		if (dcmp(Cross(B-A,C-D))){
			P ip=GetLineIntersection(A,B-A,C,D-C);
			if (OnSegment(ip,C,D)) newpoly.pb(ip);
		}
	}
	return newpoly;
}
double PolygonArea(Polygon &p){
	double area=0;
	int n=p.size();
	For(i,n-2) area+=Cross(p[i]-p[0],p[i+1]-p[0]);
	return area/2;
}
bool Onleft(Line L,P p){
	return Cross(L.v,p-L.p);
}
P GetIntersection(Line a,Line b){
	V u=a.p-b.p;
	double t=Cross(b.v,u)/Cross(a.v,b.v);
	return a.p+a.v*t;
}
int HalfplaneIntersection(Line *L,int n,P* poly){
	sort(L,L+n);
	int fi,la;
	P *p=new P[n];
	Line *q=new Line[n];
	q[fi=la=0]=L[0];
	For(i,n-1){
		while(fi<la&&!Onleft(L[i],p[la-1])) la--;
		while(fi<la&&!Onleft(L[i],p[fi])) fi++;
		q[++la]=L[i];
		if (fabs(Cross(q[la].v,q[la-1].v))<eps){
			la--;
			if (Onleft(q[la],L[i].p)) q[la]=L[i];
		}
		if (fi<la) p[la-1]=GetIntersection(q[la-1],q[la]);
	}
	while(fi<la&&!Onleft(q[fi],p[la-1])) la--;
	if (la-fi<=1) return 0;
	p[la]=GetIntersection(q[la],q[fi]);
	int m=0;
	Fork(i,fi,la) poly[m++]=p[i];
	return m;
}
P dat[4]={P(-0.5,-0.5),P(0.5,-0.5),P(0.5,0.5),P(-0.5,0.5)};
int c[4]={1,3,6,4};
int main()
{
	freopen("in.txt","r",stdin);
	freopen("out.txt","w",stdout);
	while(1)
	{
		Polygon poly;
		for (int i=0;i<4;i++){
			P A;
			if (scanf("%lf%lf",&A.x,&A.y)!=2) return 0;
			poly.pb(A);
		}
		Polygon A[4];
		for (int i=0;i<4;i++){
			A[i]=CutPolygon(poly,dat[(i+1)%4],dat[i]);
		}
		double p=0;
		for (int i=0;i<4;i++){
			p+=PolygonArea(A[i])*5/124*c[i];
		}
		p+=5.0*(5*5*4)/124;
		printf("%.10lf\n",p);
	}
	return 0;
}

//博航大腿神板
#include<iostream>
#include<cmath>
#include<algorithm>
#include<vector>
using namespace std;
#define EPS 1e-10
#define INF 1e10
#define PI 3.14159265358979323846
inline bool EQUAL(double t1, double t2){
	return t1 - t2 < EPS && t1 - t2 > -EPS;
}
inline bool LESS(double t1, double t2){
	return t1 <= t2 - EPS;
}
inline bool LESS_EQUAL(double t1, double t2){
	return t1 < t2 + EPS;
}
inline int SGN(double t){
	return LESS(t, 0) ? -1 : LESS(0, t) ? 1 : 0;
}
class Point
{
public:
	double x, y;
	Point(){}
	Point(double x, double y) :x(x), y(y){}

	bool operator == (const Point& p)const{
		return EQUAL(x, p.x) && EQUAL(y, p.y);
	}
	bool operator < (const Point& p)const{
		return LESS_EQUAL(x, p.x) && (LESS(x, p.x) || LESS(y, p.y));
	}
	Point operator + (const Point& p)const{
		return Point(x + p.x, y + p.y);
	}
	Point operator - (const Point& p)const{
		return Point(x - p.x, y - p.y);
	}
	double operator * (const Point& p)const{
		return x*p.y - y*p.x;
	}
	Point operator * (double value)const{
		return Point(x*value, y*value);
	}
	Point operator / (double value)const{
		return Point(x / value, y / value);
	}
	double dot(const Point& p)const{
		return x*p.x + y*p.y;
	}
	double r2()const{ return x*x + y*y; }
	double r()const{ return hypot(x, y); }
	double dis2(const Point& p)const{
		return (*this - p).r2();
	}
	double dis(const Point& p)const{
		return (*this - p).r();
	}

	bool onLine(const Point& p1, const Point& p2)const{
		return EQUAL((*this - p1)*(*this - p2), 0);
	}
	bool onLineSeg(const Point& p1, const Point& p2)const{
		//include extream points
		return onLine(p1, p2) && inRect(p1, p2);
	}
	double lineRelation(const Point& p1, const Point& p2)const{
		Point t = p2 - p1;
		return t.dot(*this - p1) / t.r2();
		//ret 0, *this=p1; ret 1,*this=p2;
		//ret (0,1), *this is interior to p1p2
	}
	Point footPoint(const Point& p1, const Point& p2)const{
		double r = lineRelation(p1, p2);
		return p1 + (p2 - p1)*r;
	}
	double lineDis(const Point& p1, const Point& p2)const{
		return abs((p1 - *this)*(p2 - *this)) / p1.dis(p2);
	}
	double lineSegDis(const Point& p1, const Point& p2, Point& ret)const;
	double lineSegArrayDis(const Point* p, int lineNum, Point& ret)const;
	bool lineSegArrayDisCmp(const Point* p, int lineNum, double value)const;
	Point mirror(Point& p1, Point& p2){
		Point foot = footPoint(p1, p2);
		return foot * 2 - *this;
	}

	Point rotate(double angle)const{
		Point f(sin(angle), cos(angle));
		return Point(*this * f, dot(f));
	}
	Point rotate90()const{
		return Point(-y, x);
	}
	double cosAngle(const Point& p1, const Point& p2)const{
		Point t1 = *this - p1, t2 = *this - p2;
		return t1.dot(t2) / sqrt(t1.r2()*t2.r2());
	}
	double tanAngle(const Point& o = Point(0, 0))const{
		if (EQUAL(x, o.x))return y - o.y >= 0 ? INF : -INF;
		return (y - o.y) / (x - o.x);
	}
	double angle(const Point& p1, const Point& p2)const{
		return acos(cosAngle(p1, p2));
	}
	double angle(const Point& o = Point(0, 0))const{
		return atan2(y - o.y, x - o.x);
	}
	//left return 1, right return -1, on line return 0.
	int direction(const Point& p1, const Point& p2)const{
		return SGN(x*(p1.y - p2.y) + p1.x*(p2.y - y) + p2.x*(y - p1.y));
	}
	
	bool inRect(const Point& p1, const Point& p2)const{
		return LESS_EQUAL((p1.x - x)*(p2.x - x), 0) && LESS_EQUAL((p1.y - y)*(p2.y - y), 0);
	}
	int inPolygon(const Point* p, int n)const;
	int inConvex(const Point* p, int n)const;
	int inCircle(const Point& o, double r)const{
		double dist = dis2(o);
		return SGN(r*r - dist);
	}
	void pointcut(const Point& o, double r, Point& ret1, Point& ret2)const;
	Point nearnestPoint(const Point& o, double r)const;
};
double Point::lineSegDis(const Point& p1, const Point& p2, Point& ret)const
{
	double r = lineRelation(p1, p2);
	if (LESS_EQUAL(r, 0))ret = p1;
	else if (LESS_EQUAL(1, r))ret = p2;
	else ret = footPoint(p1, p2);
	return dis(ret);
}
//input lineNum+1 points
double Point::lineSegArrayDis(const Point* p, int lineNum, Point& ret)const
{
	Point tp;
	double td, mind = INF;
	for (int i = 0; i < lineNum; i++){
		td = lineSegDis(p[i], p[i + 1], tp);
		if (LESS(td, mind)){
			mind = td; ret = tp;
		}
	}
	return mind;
}
//input lineNum+1 points
bool Point::lineSegArrayDisCmp(const Point* p, int lineNum, double value)const
{
	Point tp;
	double td;
	int flag = 1;
	for (int i = 0; i < lineNum; i++){
		td = lineSegDis(p[i], p[i + 1], tp);
		if (LESS_EQUAL(td, value))
			return true;
	}
	return false;
}

//donnot include extream points, and donnot include coincidence.
inline bool lineSegLineSegIntersect(const Point& p1, const Point& p2, const Point& q1, const Point& q2){
	Point pq1 = p1 - q1, p12 = p2 - p1, q12 = q2 - q1;
	return SGN(pq1*q12)*SGN((p2 - q1)*q12) < 0 && SGN(pq1*p12)*SGN((p1 - q2)*p12) < 0;
}
//include extream points and coincidence.
inline bool lineSegLineSegIntersect2(const Point& p1, const Point& p2, const Point& q1, const Point& q2){
	if (!(LESS_EQUAL(min(q1.x, q2.x), max(p1.x, p2.x)) && LESS_EQUAL(min(p1.x, p2.x), max(q1.x, q2.x))
		&& LESS_EQUAL(min(q1.y, q2.y), max(p1.y, p2.y)) && LESS_EQUAL(min(p1.y, p2.y), max(q1.y, q2.y))))
		return false;
	Point pq1 = p1 - q1, p12 = p2 - p1, q12 = q2 - q1;
	return SGN(pq1*q12)*SGN((p2 - q1)*q12) <= 0 && SGN(pq1*p12)*SGN((p1 - q2)*p12) <= 0;
}
//donot include extream points, and donot include coincidence.
inline bool lineLineSegIntersect(const Point& l1, const Point& l2, const Point& p1, const Point& p2){
	Point line = l2 - l1;
	return SGN((p1 - l1)*line)*SGN((p2 - l1)*line) < 0;
}
//donnot include coincidence.
inline bool lineLineIntersect(const Point& p1, const Point& p2, const Point& q1, const Point& q2){
	return !EQUAL((p2 - p1)*(q2 - q1), 0);
}
inline Point lineLineIntersectPoint(const Point& p1, const Point& p2, const Point& q1, const Point& q2){
	Point q12 = q2 - q1;
	double k = (p2 - p1)*q12;
	if (EQUAL(k, 0))return Point(INF*INF, INF*INF);
	double r = ((q1 - p1)*q12) / k;
	return p1 + (p2 - p1) * r;
}

Point circumcenter(const Point& p1, const Point& p2, const Point& p3)
{
	Point t1 = (p1 + p2)*0.5, t2, t3 = (p2 + p3)*0.5, t4;
	t2 = t1 + (p1 - p2).rotate90();
	t4 = t3 + (p2 - p3).rotate90();
	return lineLineIntersectPoint(t1, t2, t3, t4);
}
Point incenter(const Point& p1, const Point& p2, const Point& p3)
{
	double r12 = p1.dis(p2), r23 = p2.dis(p3), r31 = p3.dis(p1);
	Point t1 = (p2*r31 + p3*r12) / (r12 + r31), t2 = (p1*r23 + p3*r12) / (r12 + r23);
	return lineLineIntersectPoint(p1, t1, p2, t2);
}
Point prepencenter(const Point& p1, const Point& p2, const Point& p3)
{
	Point t1 = p1 + (p2 - p3).rotate90();
	Point t2 = p2 + (p1 - p3).rotate90();
	return lineLineIntersectPoint(p1, t1, p2, t2);
}
inline Point barycenter(const Point& p1, const Point& p2, const Point& p3){
	return (p1 + p2 + p3) / 3;
}
inline double apothem(const Point& p1, const Point& p2, const Point& p3){
	Point p12 = p2 - p1, p13 = p3 - p1, p23 = p3 - p2;
	return abs(p12*p23) / (p12.r() + p13.r() + p23.r());
}
inline double circumradius(const Point& p1, const Point& p2, const Point& p3){
	Point p12 = p2 - p1, p13 = p3 - p1, p23 = p3 - p2;
	return sqrt(p12.r2()*p23.r2()*p13.r2()) / (2 * abs(p12*p23));
}

int getPolygonDirection(const Point* p, int n)
{
	int index = 0;
	for (int i = 1; i < n; i++){
		if (p[i] < p[index])index = i;
	}
	return p[index].direction(p[index + 1 < n ? index + 1 : 0], p[index - 1 >= 0 ? index - 1 : n - 1]);
}
bool checkConvex(const Point* p, int n)
{
	int direction = p[0].direction(p[n - 1], p[1]);
	if (direction == 0)return false;
	if (p[n - 1].direction(p[n - 2], p[0]) != direction)return false;
	for (int i = n - 2; i > 0; i--){
		if (p[i].direction(p[i - 1], p[i + 1]) != direction)
			return false;
	}
	return true;
}
bool checkConvex(const Point* p, int n, bool *ret)
{
	bool retValue = true;
	int direction = getPolygonDirection(p, n);
	if (!(ret[n - 1] = p[n - 1].direction(p[0], p[n - 2]) == direction))
		retValue = false;
	if (!(ret[0] = p[0].direction(p[1], p[n - 1]) == direction))
		retValue = false;
	for (int i = n - 2; i > 0; i--){
		if (!(ret[i] = p[i].direction(p[i + 1], p[i - 1]) == direction))
			retValue = false;
	}
	return retValue;
}
double polygonArea(const Point* p, int n)
{
	double area = 0;
	for (int i = n - 2; i > 0; i--)
		area += p[i].y *(p[i - 1].x - p[i + 1].x);
	area += p[0].y*(p[n - 1].x - p[1].x);
	area += p[n - 1].y*(p[n - 2].x - p[0].x);
	return area / 2;
}
int Point::inPolygon(const Point* p, int n)const
{
	int i, j = n - 1, odd = -1;
	for (i = 0; i < n; j = i++){
		if (LESS(p[i].y, y) != LESS(p[j].y, y)){
			double tx = (y - p[j].y) / (p[i].y - p[j].y)*(p[i].x - p[j].x) + p[j].x;
			if (LESS_EQUAL(tx, x)){
				if (LESS(tx, x))odd = -odd;
				else return 0;
			}
		}
		else if (onLineSeg(p[i], p[j]))return 0;
	}
	return odd;
}
int Point::inConvex(const Point* p, int n)const
{
	int _direction = p[1].direction(p[2], p[0]);
	if (direction(p[0], p[1]) != _direction){
		if (onLineSeg(p[0], p[1]))return 0;
		return -1;
	}
	if (direction(p[n - 1], p[0]) != _direction){
		if (onLineSeg(p[n - 1], p[0]))return 0;
		return -1;
	}
	int left = 2, right = n - 1;
	while (left < right){
		int mid = (left + right) >> 1;
		if (direction(p[0], p[mid]) == _direction)left = mid + 1;
		else right = mid;
	}
	int ret = direction(p[left-1],p[left]);
	return ret == _direction ? 1 : ret == 0 ? 0 : -1;
}
Point lineConvexIntersectPointInternal(const Point& p1, const Point& p2, const Point* p, int n, int start, int end)
{
	Point p12 = p2 - p1;
	if (end < start)end += n;
	double value = SGN((p[start] - p1)*p12);
	while (start + 1 < end){
		int mid = (start + end) / 2;
		Point cur = p[mid < n ? mid : mid - n];
		double t = (cur - p1)*p12*value;
		if (LESS(0, t))start = mid;
		else if (LESS(t, 0))end = mid;
		else return cur;
	}
	if (start >= n)start -= n;
	return lineLineIntersectPoint(p1, p2, p[start], p[start + 1]);
}
int lineConvexIntersectPoint(const Point& p1, const Point& p2, const Point* p, int n, Point& ret1, Point& ret2)
{
	Point p12 = p2 - p1;
	int pos = 0, step = n * 2 / 3;
	double d = (p[pos] - p1)*p12;
	int zero = -1, pos2 = -1;
	while (step > 1){
		step=(step + 1) / 2;
		int i = pos + step, k = pos - step;
		if (i >= n)i -= n;
		if (k < 0)k += n;
		double di = (p[i] - p1)*p12, dk = (p[k] - p1)*p12;
		if (SGN(di)*SGN(d) < 0){ pos2 = i; break; }
		if (SGN(dk)*SGN(d) < 0){ pos2 = k; break; }
		if (abs(di) < abs(d)){ d = di; pos = i; }
		if (abs(dk) < abs(d)){ d = dk; pos = k; }
		if (EQUAL(d, 0)){ zero = pos; break; }
	}
	if (zero != -1){
		ret1 = p[zero];
		int left = zero - 1 >= 0 ? zero - 1 : n - 1;
		int right = zero + 1 < n ? zero + 1 : 0;
		double dl = (p[left] - p1)*p12, dr = (p[right] - p1)*p12;
		if (EQUAL(dl, 0)){ ret2 = p[left]; return 3; }
		else if (EQUAL(dr, 0)){ ret2 = p[right]; return 3; }
		else if (dl*dr < 0)return 1;
		else{ pos = left; pos2 = right; }
	}
	if (pos2 == -1)return 0;
	ret1 = lineConvexIntersectPointInternal(p1, p2, p, n, pos, pos2);
	ret2 = lineConvexIntersectPointInternal(p1, p2, p, n, pos2, pos);
	return 2;
}

bool lineSegInPolygon(const Point& p1, const Point& p2, const Point* p, int n)
{
	bool flag = false;
	Point minPoint;
	switch (p1.inPolygon(p, n)){
	case -1:return false;
	case 0:flag = true;
	}
	switch (p2.inPolygon(p, n)){
	case -1:return false;
	case 1:flag = false;
	}
	if (flag)minPoint = max(p1, p2);
	for (int i = 0, j = n - 1; i < n; j = i++){
		if (p[i].onLineSeg(p1, p2) && !(p[i] == p1 || p[i] == p2)){
			if (p[i > 0 ? i - 1 : n - 1].direction(p1, p2) * p[i + 1 < n ? i + 1 : 0].direction(p1, p2) < 0)
				return false;
			if (flag && p[i] < minPoint)minPoint = p[i];
		}
		else if (lineSegLineSegIntersect(p[i], p[j], p1, p2))
			return false;
	}
	if (flag){
		const Point& t = min(p1, p2);
		Point mid = (t + minPoint)*0.5;
		if (mid.inPolygon(p, n) == -1)return false;
	}
	return true;
}
Point gravityCenter(const Point* p, int n)
{
	if (n < 3){
		if (n == 1)return p[0];
		else return (p[0] + p[1])*0.5;
	}
	double area = 0;
	Point ret(0, 0);
	for (int i = 0, j = n - 1; i < n; j = i++){
		double t = p[i] * p[j];
		area += t;
		ret.x += (p[i].x + p[j].x)*t;
		ret.y += (p[i].y + p[j].y)*t;
	}
	return ret / (3 * area);
}
//ret[n] must be available to visit.
int convexHullSorted(const Point* p, int n, Point* ret)
{
	int j = 0;
	for (int i = 0; i < n; i++){
		while (j >= 2 && p[i].direction(ret[j - 2], ret[j - 1]) != 1)j--;
		ret[j++] = p[i];
	}
	int mid = j + 1;
	for (int i = n - 2; i >= 0; i--){
		while (j >= mid && p[i].direction(ret[j - 2], ret[j - 1]) != 1)j--;
		ret[j++] = p[i];
	}
	return j - 1;
}
void convexHullSorted(const Point* p, int n, Point* up, int& retUp, Point* down, int& retDown)
{
	retUp = retDown = 0;
	for (int i = 0; i < n; i++){
		while (retUp >= 2 && p[i].direction(up[retUp - 2], up[retUp - 1]) != -1)retUp--;
		while (retDown >= 2 && p[i].direction(down[retDown - 2], down[retDown - 1]) != 1)retDown--;
		up[retUp++] = p[i];
		down[retDown++] = p[i];
	}
}
int halfPlainIntersectInternal(vector<pair<double, const Point*>>& v, int n, Point* ret)
{
	for (int i = 0; i < n; i++)
		v[i].first = v[i].second[1].angle(v[i].second[0]);
	sort(v.begin(), v.end());
	vector<const Point*> line(n);
	vector<Point> point(n);
	int first = 0, last = 0;
	line[0] = v[0].second;
	for (unsigned int i = 1; i < v.size(); i++){
		while (first < last && point[last - 1].direction(v[i].second[0], v[i].second[1]) == -1) last--;
		while (first < last && point[first].direction(v[i].second[0], v[i].second[1]) == -1) first++;
		line[++last] = v[i].second;
		if (!lineLineIntersect(line[last - 1][0], line[last - 1][1], line[last][0], line[last][1])){
			last--;
			if (v[i].second[0].direction(line[last][0], line[last][1]) == 1)line[last] = v[i].second;
		}
		if (first<last)
			point[last - 1] = lineLineIntersectPoint(line[last - 1][0], line[last - 1][1], line[last][0], line[last][1]);
	}
	while (first < last && point[last - 1].direction(line[first][0], line[first][1]) == -1) last--;
	if (last - first <= 1) return 0;
	point[last] = lineLineIntersectPoint(line[first][0], line[first][1], line[last][0], line[last][1]);
	int num = unique(&*point.begin() + first, &*point.begin() + last + 1) - &point[first];
	while (num>1 && point[first] == point[first + num - 1])num--;
	memcpy(ret, &point[first], sizeof(Point)*num);
	return num;
}
int halfPlainIntersect(const Point(*p)[2], int n, Point* ret)
{
	vector<pair<double, const Point*>> v(n + 4);
	Point ext[4][2] = { { { -INF, -INF }, { INF, -INF } }, { { INF, -INF }, { INF, INF } },
	{ { INF, INF }, { -INF, INF } }, { { -INF, INF }, { -INF, -INF } } };
	for (int i = 0; i < 4; i++)
		v[i].second = ext[i];
	for (int i = 0; i < n; i++)
		v[i + 4].second = p[i];
	return halfPlainIntersectInternal(v, n + 4, ret);
}
int polygonKernel(const Point* p, int n, Point* ret)
{
	vector<pair<double, const Point*>> v;
	Point ext[2] = { p[n - 1], p[0] };
	v[0].second = ext;
	for (int i = 1; i < n; i++)
		v[i].second = &p[i - 1];
	return halfPlainIntersectInternal(v, n, ret);
}

struct NearestPointsStruct{
	Point p1, p2;
	double d2;
	vector<Point> v;
};
inline bool nearestPointsCmp(const Point& p1, const Point& p2){
	return LESS_EQUAL(p1.y, p2.y) && (LESS(p1.y, p2.y) || LESS(p1.x, p2.x));
}
void nearestPointsInternal(const Point* p, int left, int right, NearestPointsStruct& s)
{
	if (right - left < 8){
		for (int i = left; i < right; i++){
			for (int j = i + 1; i < right; j++){
				double td2 = p[j].dis2(p[i]);
				if (td2 < s.d2){
					s.d2 = td2;
					s.p1 = p[i]; s.p2 = p[j];
				}
			}
		}
		return;
	}
	int mid = (left + right) >> 1;
	nearestPointsInternal(p, left, mid, s);
	nearestPointsInternal(p, mid, right, s);
	s.v.clear();
	double l = (p[mid - 1].x + p[mid].x) / 2;
	for (int i = mid - 1; i >= left && (p[i].x - l)*(p[i].x - l) < s.d2; i++)
		s.v.push_back(p[i]);
	for (int i = mid; i<right && (p[i].x - l)*(p[i].x - l) < s.d2; i++)
		s.v.push_back(p[i]);
	sort(s.v.begin(), s.v.end(), nearestPointsCmp);
	for (unsigned int i = 0; i < s.v.size(); i++){
		for (unsigned int j = i + 1; j < s.v.size() && (p[j].y - p[i].y)*(p[j].y - p[i].y) < s.d2; j++){
			double td2 = p[j].dis2(p[i]);
			if (td2 < s.d2){
				s.d2 = td2;
				s.p1 = p[i]; s.p2 = p[j];
			}
		}
	}
}
double nearestPointsSorted(const Point* p, int n, Point& ret1, Point& ret2)
{
	NearestPointsStruct s;
	s.d2 = INF;
	s.v.reserve(n);
	nearestPointsInternal(p, 0, n, s);
	ret1 = s.p1; ret2 = s.p2;
	return sqrt(s.d2);
}
double farthestPointsConvex(const Point* p, int n, Point& ret1, Point& ret2)
{
	double d2 = 0;
	for (int i = n - 1, j = n - 2; i >0; i--){
		while (1){
			double td2 = p[i].dis2(p[j]);
			if (td2 > d2){
				d2 = td2;
				ret1 = p[i]; ret2 = p[j];
			}
			if (!j)break;
			j--;
		}
	}
	return sqrt(d2);
}
double farthestPointsSorted(const Point* p, int n, Point& ret1, Point& ret2)
{
	vector<Point> v;
	v.reserve(n);
	//convexHullSorted(p, n, &*v.begin());
	return farthestPointsConvex(&*v.begin(), v.size(), ret1, ret2);
}

int circleLineRelation(const Point& o, double r, const Point& p1, const Point& p2)
{
	double d = o.lineDis(p1, p2);
	if (LESS(d, r))return 1;
	if (LESS(r, d))return 3;
	return 2;
}
int circleCircleRelation(const Point& o1, double r1, const Point& o2, double r2)
{
	double r = o1.dis(o2);
	if (LESS(r1 + r2, r))return 4;
	if (!LESS_EQUAL(r1 + r2, r))return 3;
	double sub = abs(r1 - r2);
	if (LESS(sub, r))return 2;
	if (!LESS_EQUAL(sub, r))return 1;
	return 0;
}
bool circleLineSegIntersect(const Point& o, double r, const Point& p1, const Point& p2)
//include extream points.
{
	int t1 = p1.inCircle(o, r), t2 = p2.inCircle(o, r);
	if (t1 >= 0 || t2 >= 0)
		return t1 != 1 || t2 != 1;
	double t = o.lineRelation(p1, p2);
	if (t >= 1 || t <= 0)return false;
	Point foot = p1 + (p2 - p1)*r;
	return foot.inCircle(o, r) >= 0;
}
void circleLineIntersect(const Point& o, double r, const Point& p1, const Point& p2, Point& ret1, Point& ret2)
{
	Point foot = o.footPoint(p1, p2);
	double t = sqrt((r*r - o.dis2(foot)) / p1.dis2(p2));;
	ret1 = foot + (p2 - p1)*t;
	ret2 = foot * 2 - ret1;
}
void circleCircleIntersect(const Point& o1, double r1, const Point& o2, double r2, Point& ret1, Point& ret2)
{
	double d2 = o1.dis2(o2);
	double t1 = (r1*r1 - r2*r2) / (2 * d2) + 0.5;
	double t2 = sqrt(r1*r1 / d2 - t1*t1);
	Point foot = o1 + (o2 - o1)*t1;
	ret1 = foot + (o2 - o1).rotate90()*t2;
	ret2 = foot * 2 - ret1;
}
void Point::pointcut(const Point& o, double r, Point& ret1, Point& ret2)const
{
	double t1 = r*r / dis2(o);
	Point foot = o + (o - *this)*t1;
	double t2 = sqrt(t1 - t1*t1);
	ret1 = foot + (*this - o).rotate90()*t2;
	ret2 = foot * 2 - ret1;
}
Point Point::nearnestPoint(const Point& o, double r)const
{
	Point p = *this - o;
	double d = p.r();
	if (EQUAL(d, 0))return o;
	return o + p*(r / d);
}
//Upset the order before using this function.
double minCoveringCircle(const Point* p, int n, Point& ret)
{
	if (n == 1){ ret = p[0]; return 0; }
	double r2 = p[0].dis2(p[1]);
	ret = (p[0] + p[1]) * 0.5;
	for (int i = 2; i < n; i++){
		if (LESS(r2, ret.dis2(p[i]))){
			ret = (p[0] + p[i]) * 0.5;
			r2 = p[0].dis2(p[i]);
			for (int j = 1; j < i; j++){
				if (LESS(r2, ret.dis2(p[j]))){
					ret = (p[i] + p[j]) * 0.5;
					r2 = p[i].dis2(p[j]);
					for (int k = 0; k < j; k++){
						if (LESS(r2, ret.dis2(p[k]))){
							ret = circumcenter(p[i], p[j], p[k]);
							r2 = ret.dis2(p[k]);
						}
					}
				}
			}
		}
	}
	return sqrt(r2);
}
int unitCoveringCircle(const Point* p, int n, double r)
{
	int ret = 0;
	vector<pair<double, bool>> v;
	v.reserve(2 * n);
	double t = r*r * 4;
	for (int i = 0; i < n; i++){
		v.clear();
		int value = 0;
		for (int j = 0; j < n; j++){
			if (LESS_EQUAL(p[i].dis2(p[j]), t) && i != j){
				double a = p[j].angle(p[i]);
				double b = acos(p[i].dis(p[j]) / r / 2);
				double t1 = a - b, t2 = a + b;
				if (t1 < -PI / 2){
					if (t2 < -PI / 2){
						a += 2 * PI;
						b += 2 * PI;
					}
					else value++;
				}
				v.push_back(make_pair(t1, true));
				v.push_back(make_pair(t2, false));
			}
		}
		sort(v.begin(), v.end());
		if (value > ret)ret = value;
		for (unsigned int j = 0; j < v.size(); j++){
			if (v[j].second){
				value++;
				if (value > ret)ret = value;
			}
			else value--;
		}
	}
	return ret;
}

```

