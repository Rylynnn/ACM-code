#include<iostream>
#include<cstdio>
#include<queue>
#define maxn 1005					//µãÊý
#define maxm 80005					//±ßÊý
#define INF 0x3f3f3f3f
#define rever(x) (mem+((x-mem)^1))
using namespace std;
struct edge
{
	int s,t,v,c;
	edge* next;
}mem[maxm],*head[maxn],*prev[maxn];
queue<int> q;
int cnt=-1,n;
int dis[maxn];
int S,T;
void add_edge(int s,int t,int v,int c)
{
	mem[++cnt].s=s;mem[cnt].t=t;mem[cnt].v=v;mem[cnt].c=c;mem[cnt].next=head[s];head[s]=mem+cnt;
	mem[++cnt].s=t;mem[cnt].t=s;mem[cnt].v=0;mem[cnt].c=-c;mem[cnt].next=head[t];head[t]=mem+cnt;
}
bool bfs()
{
	for (int i=0;i<=n;i++) dis[i]=INF;
	q.push(S);dis[S]=0;
	while(!q.empty())
	{
		for (edge *it=head[q.front()];it;it=it->next)
		if (it->v&&dis[q.front()]+it->c<dis[it->t])
		{
			dis[it->t]=dis[q.front()]+it->c;
			prev[it->t]=it;
			q.push(it->t);
		}
		q.pop();
	}
	return (dis[T]!=INF);
}
int cost=0;
int dinic()
{
	int flow=0;
	while(bfs())
	{
		int augflow=INF,tmpcost=0;
		for (edge* it=prev[T];it;it=prev[it->s])
		{
			augflow=min(augflow,it->v);
			tmpcost+=it->c;
		}
		for (edge* it=prev[T];it;it=prev[it->s])
		{
			it->v-=augflow;
			rever(it)->v+=augflow;
		}
		flow+=augflow;cost+=augflow*tmpcost;
	}
	return flow;
}
int N,M,A,B,C;
int main()
{
	scanf("%d%d",&N,&M);
	S=0;T=N+1;n=T;
	add_edge(S,1,2,0);add_edge(N,T,2,0);
	for (int i=1;i<=M;i++)
	{
		scanf("%d%d%d",&A,&B,&C);
		add_edge(A,B,1,C);
		add_edge(B,A,1,C);
	}
	dinic();
	printf("%d\n",cost);
	return 0;
}

