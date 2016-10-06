int nxt[maxn];
char origin_string[maxn];
char target_string[maxn];
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
int main()
{
	int _;
	scanf("%d",&_);
	while(_--)
	{
		scanf("%s",target_string);
		scanf("%s",origin_string);
		get_nxt();
		printf("%d\n",kmp());
	}
	return 0;
}
