#include<iostream>
#include<cstdio>
#include<queue>
#define INF 0x3f3f3f3f
using namespace std;
struct edge
{
	int n,v;
	edge* next;
}mem[2000005],*head[2000],*store[505];
queue<int> q;
int cnt,n;
int dis[2000];
char ch;
void add_edge(int s,int t,int v)
{
	mem[cnt].n=t;mem[cnt].v=v;
	mem[cnt].next=head[s];head[s]=mem+cnt;
	cnt++;
	mem[cnt].n=s;mem[cnt].v=0;
	mem[cnt].next=head[t];head[t]=mem+cnt;
	cnt++;
}
bool bfs(int x)
{
	for (int i=0;i<=n;i++) dis[i]=INF;
	q.push(x);dis[x]=0;
	while(!q.empty())
	{
		for (edge *it=head[q.front()];it;it=it->next)
		if (it->v&&dis[it->n]==INF)
		{
			dis[it->n]=dis[q.front()]+1;
			q.push(it->n);
		}
		q.pop();
	}
	return (dis[n]!=INF);
}
int dinic(int x,int f)
{
	if (x==n) return f;
	for (edge* it=head[x];it;it=it->next)
	if (dis[it->n]==dis[x]+1&&it->v)
	{
		int nowf=dinic(it->n,min(f,it->v));
		if (nowf)
		{
			it->v-=nowf;
			mem[(it-mem)^1].v+=nowf;
			return nowf;
		}
	}
	return 0;
}
int DINIC(int x)
{
	int ths=0;
	while(bfs(x))
	{
		while(1)
		{
			int tmp=dinic(x,INF);
			if (!tmp) break;
			ths+=t;
		}
	}
	return ths;
}
//flow=DINIC(source);
