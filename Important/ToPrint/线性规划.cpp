这是一道模板题。

本题中你需要求解一个标准型线性规划：

有nn个实数变量x1,x2,⋯,xnx1,x2,⋯,xn和mm条约束，其中第ii条约束形如∑nj=1aijxj≤bi∑j=1naijxj≤bi。

此外这nn个变量需要满足非负性限制，即xj≥0xj≥0。

在满足上述所有条件的情况下，你需要指定每个变量xjxj的取值，使得目标函数F=∑nj=1cjxjF=∑j=1ncjxj的值最大。

输入格式
第一行三个正整数 n,m,tn,m,t。其中t∈{0,1}t∈{0,1}。

第二行有nn个整数c1,c2,⋯,cnc1,c2,⋯,cn，整数间均用一个空格分隔。

接下来mm行，每行代表一条约束，其中第ii行有n+1n+1个整数ai1,ai2,⋯,ain,biai1,ai2,⋯,ain,bi，整数间均用一个空格分隔。

输出格式
如果不存在满足所有约束的解，仅输出一行"Infeasible"。

如果对于任意的MM，都存在一组解使得目标函数的值大于MM，仅输出一行"Unbounded"。

否则，第一行输出一个实数，表示目标函数的最大值FF。当第一行与标准答案的相对误差或绝对误差不超过10−610−6，你的答案被判为正确。

如果t=1t=1，那么你还需要输出第二行，用空格隔开的nn个非负实数，表示此时x1,x2,⋯,xnx1,x2,⋯,xn的取值，如有多组方案请任意输出其中一个。

判断第二行是否合法时，我们首先检验F−∑nj=1cjxjF−∑j=1ncjxj是否为00，再对于所有ii，检验min{0,bi−∑nj=1aijxj}min{0,bi−∑j=1naijxj}是否为00。检验时我们会将其中大于00的项和不大于00的项的绝对值分别相加得到S+S+和S−S−，如果S+S+和S−S−的相对误差或绝对误差不超过10−610−6，则判为正确。

如果t=0t=0，或者出现Infeasible或Unbounded时，不需要输出第二行。

//Love LBL
//Author:zhouzixuan
//单纯形 
//uoj179
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cmath>
using namespace std;
const int N=30;
const double zero=1e-8;
const double INF=1e15;
int n,m,type; double a[N][N],ans[N];
inline long long getint()
{
    long long x=0; char c=getchar(); bool flag=false;
    while ((c!='-')&&((c<'0')||(c>'9'))) c=getchar();
    if (c=='-') flag=true,c=getchar();
    while ((c>='0')&&(c<='9')) x=x*10+(long long)(c-'0'),c=getchar();
    if (flag) return -x; else return x;
}
void init()
{
	n=getint(); m=getint(); type=getint();
	for (int i=1; i<=n; i++) a[0][i]=getint();
	for (int i=1; i<=m; i++)
	{
		for (int j=1; j<=n; j++) a[i][j]=getint();
		a[i][0]=getint();
	}
}
namespace Simplex
{
	int tot,b[N],idx[N],idy[N];
	inline int dcmp(double x)
	{
		if (x>zero) return 1;
		if (x<-zero) return -1;
		return 0;
	}
	void pivot(int x,int y)
	{
		swap(idy[x],idx[y]);
		tot=0; double tmp=a[x][y]; a[x][y]=1/a[x][y];
		for (int i=0; i<=n; i++) if (y!=i) a[x][i]/=tmp;
		for (int i=0; i<=n; i++) if (dcmp(a[x][i])!=0) b[++tot]=i;
		for (int i=0; i<=m; i++)
		{
			if ((x==i)||(dcmp(a[i][y])==0)) continue;
			for (int j=1; j<=tot; j++) if (b[j]!=y) a[i][b[j]]-=a[x][b[j]]*a[i][y];
			a[i][y]=-a[i][y]/tmp;
		}
	}
	bool init()
	{
		for (int i=1; i<=n; i++) idx[i]=i;
		for (int i=1; i<=m; i++) idy[i]=i+n;
		while (true)
		{
			int x=0,y=0;
			for (int i=1; i<=m; i++) if ((dcmp(a[i][0])<0)&&((x==0)||(rand()&1))) x=i; if (x==0) break;
			for (int i=1; i<=n; i++) if ((dcmp(a[x][i])<0)&&((y==0)||(rand()&1))) y=i;
			if (y==0) return false; pivot(x,y);
		}
		return true;
	}
	void solve()
	{
		if (!init()) {printf("Infeasible\n"); return;}
		while (true)
		{
			int x=0,y=0; double mint=INF;
			for (int i=1; i<=n; i++) if (dcmp(a[0][i])>0) {x=i; break;} if (x==0) break;
			for (int i=1; i<=m; i++) if ((dcmp(a[i][x])>0)&&(a[i][0]/a[i][x]<mint)) mint=a[i][0]/a[i][x],y=i;
			if (y==0) {printf("Unbounded\n"); return;} pivot(y,x);
		}
		printf("%.8lf\n",-a[0][0]); if (type==0) return;
		for (int i=1; i<=m; i++) if (idy[i]<=n) ans[idy[i]]=a[i][0];
		for (int i=1; i<=n; i++) printf("%.8lf ",ans[i]); printf("\n");
	}
}
int main()
{
	srand(2333333);
	init();
	Simplex::solve();
	return 0;
}