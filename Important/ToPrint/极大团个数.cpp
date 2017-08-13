//邻接矩阵输入
//bitset优化
//ans为极大团个数
//本题点数不超过40

/*
 * 原理
 * 这个算法主要是构造了三个集合
R集合记录的是当前极大团中已经加入的点
P集合记录的是可能还能加入的点（也就是说可能与R集合中所有点都有边存在的点）
X集合记录的是已经完成极大团计数的点（作用是判重）
P ∪ X是所有可能与R集合构成极大团的点集（虽然我们已经知道X中的点不可能在参与极大团的构成），也就是与最后一个加入R集合相连的点相连的点的一部分
初始是RX都是空集合，P中是所有点的集合，那么如何对一个个点进行操作呢
   BronKerbosch1(R, P, X):
       if P and X are both empty:
           report R as a maximal clique
       for each vertex v in P:
           BronKerbosch1(R ⋃ {v}, P ⋂ N(v), X ⋂ N(v))
           P := P \ {v}
           X := X ⋃ {v}
这是基础的Born_Kerbosch算法，对于每一个点P中的点v我们把v加入集合R，对在P集合中与点v相连的一部分中寻找下一个可能加入R集合的点，然后我们把v从P集合中移出，加入X集合代表当前状态下对包含点v的极大团已经计算完毕了。R集合为极大团的时候，必须要满足P与X都是空的，P存放的是还可能加入R集合的点，P集合为空代表没有点还能再加入到R集合当中，而X集合存放的是已经完成极大团计数的点，而且X集合中的点必然是与所有R集合中的点都有边存在的，也即X集合中点必然可以与R集合构成极大团，如果X集合不是空的的话，那么说明R集合中的极大团是在之前计算包含X集合中的点的极大团的时候已经计算过了的，故当且仅当PX都为空集合的时候R才是一个极大团，ans自增1

然后是Born_Kerbosch的一个优化
我们知道在上述的算法中必然有许多重复计算之前计算过的极大团然后调出的过程。然后我们考虑如下问题，取集合P∪X中的一个点u，要与R集合构成极大团，那么取得点必然是P \ N(u)中一个一个点 N（u）代表P中是u邻居的点，通俗的来讲取的点要么是u，要么就是与u不相邻的点，这样我们照样可以重复不漏的计算所有极大团，这个的正确性十分好想，就不证明了，也可以参照wiki中的解释，然后我们每次取从原先的遍历P中所有点变成了遍历P中u点与不与u相邻的点即P \ N(u)，这样我们可以减少许多不必要的计算，而我们要想进一步减少计算，我们就可以取u，使得u的邻居最多，也就是使我们要遍历的点进一步减少
 */

/*
 * #include<cstdio>    
#include<cstring>    
using namespace std;    
const int N=130;    
int ans,a[N][N],R[N][N],P[N][N],X[N][N];//a记录图    
bool Bron_Kerbosch(int d,int nr,int np,int nx)//d记录当前计算的是第几个点nr是R中点的个数，Np，NX依次类推    
{    
    int i,j;    
    if(np==0&&nx==0)//PX为空 输出极大团    
    {    
        ans++;    
        if(ans>1000)//超过题目限制 跳出所有循环    
            return 1;    
        return 0;    
    }    
    int u,max=0;    
    u=P[d][1];    
    for(i=1;i<=np;i++)    
    {    
        int cnt=0;    
        for(j=1;j<=np;j++)    
        {    
            if(a[P[d][i]][P[d][j]])    
                cnt++;    
        }    
        if(cnt>max)    
        {    
            max=cnt;    
            u=P[d][i];    
        }    
    }    
    for(i=1;i<=np;i++)    
    {    
        int v=P[d][i];    
        if(a[v][u]) continue;    
        for(j=1;j<=nr;j++)    
            R[d+1][j]=R[d][j];    
        R[d+1][nr+1]=v;    
        int cnt1=0;    
        for(j=1;j<=np;j++)    
            if(P[d][j]&&a[P[d][j]][v])    
                P[d+1][++cnt1]=P[d][j];    
        int cnt2=0;    
        for(j=1;j<=nx;j++)    
            if(a[X[d][j]][v])    
                X[d+1][++cnt2]=X[d][j];    
        if(Bron_Kerbosch(d+1,nr+1,cnt1,cnt2))    
            return 1;    
        P[d][i]=0;    
        X[d][++nx]=v;    
    }    
    return 0;    
}    
int main()    
{    
    int n,i,m,x,y;    
    while(scanf("%d%d",&n,&m)!=EOF)    
    {    
        memset(a,0,sizeof(a));    
        while(m--)    
        {    
            scanf("%d%d",&x,&y);    
            a[x][y]=a[y][x]=1;    
        }    
        ans=0;    
        for(i=1;i<=n;i++)    
            P[1][i]=i;    
        Bron_Kerbosch(1,0,n,0);    
        if(ans>1000)    
            printf("Too many maximal sets of friends.\n");    
        else    
            printf("%d\n",ans);    
    }    
    return 0;    
}    
 */

 //Claris
#include<cstdio>
#include<algorithm>
using namespace std;
typedef unsigned long long U;
typedef long long ll;
const int N=45;
int n,K,x,i,j,ans;bool flag;U g[N];double res;
inline int ctz(U s){return s?__builtin_ctzll(s):64;}
void BornKerbosch(U cur,U allow,U forbid){
  if(!allow&&!forbid){
    ans=max(ans,__builtin_popcountll(cur));
    return;
  }
  if(!allow)return;
  int pivot=ctz(allow|forbid);
  U z=allow&~g[pivot];
  for(int u=ctz(z);u<n;u+=ctz(z>>(u+1))+1){
    BornKerbosch(cur|(1ULL<<u),allow&g[u],forbid&g[u]);
    allow^=1ULL<<u,forbid|=1ULL<<u;
  }
}
int main(){
  scanf("%d%d",&n,&K);
  for(i=0;i<n;i++)g[i]=(1ULL<<n)-1-(1ULL<<i);
  for(i=0;i<n;i++)for(j=0;j<n;j++){
    scanf("%d",&x);
    if(!x&&i!=j)g[i]^=1ULL<<j;
  }
  BornKerbosch(0,(1ULL<<n)-1,0);
  if(ans>1)res=1.0*(ans-1)/2.0/ans;
  printf("%.10f",res*K*K);
}
