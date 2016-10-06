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
	

int main()
{
	scanf("%d\n",&__);
	for(_=1;_<=__;_++)
	{
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
	}
	return 0;
}
