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

