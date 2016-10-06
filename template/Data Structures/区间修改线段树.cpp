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

inline void read(int &x) {
    char ch=getchar();
    while(ch<'0'||ch>'9') ch=getchar();
    x=0;
    while(ch<='9'&&ch>='0'){
        x=x*10+ch-'0';
        ch=getchar();
    }
}

int main()
{
	read(n);read(m);
	for(i=1;i<=n;i++)read(c[i]);
	k=1;mt(1,n);
	for(i=1;i<=m;i++)
	{
		read(l);read(j);read(k);
		if(find(j,k,1)<l)
		{
			printf("-1\n%d",i);
			break;
		}
		add(j,k,-l,1);
	}
	if(i>m)printf("0");
	return 0;
}


