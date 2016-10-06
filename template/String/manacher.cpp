char s1[1000005],s2[2000010];
int p[2000010];
int i,j,k,l,m,n;
int mx,id;

int min(int x,int y)
{
	return x<y?x:y;
}

int main()
{
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
	return 0;
}
