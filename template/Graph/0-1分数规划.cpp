const int maxn=1005;
double a[maxn],b[maxn];
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

