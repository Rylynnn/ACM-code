#include<cstdio>
#include<cstring>

int next[200005];
int i,j,k,l,ns,nm;
char s[500005],m[200005];

int main()
{
	gets(s);
	gets(m);
	ns=strlen(s);nm=strlen(m);
	next[0]=-1;i=0;j=-1;
	while(i<nm)
	{
		if(j==-1||s[i]==s[j])
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
	return 0;
}
