#include<iostream>
#include<cstdio>
#include<queue>
#include<algorithm>
#include<cstring>
#include<cmath>
using namespace std;
#define maxn 205						//最大点数
#define maxm 205						//最大边数
#define rever(x) (mem+((x-mem)^1))
struct edge
{
	int s,t,v;
	edge* next;
}mem[maxm*2],*head[maxn];
int cnt=-1;
void add_edge(int s,int t,int v)
{
	mem[++cnt].s=s;mem[cnt].t=t;mem[cnt].v=v;mem[cnt].next=head[s];head[s]=mem+cnt;
	mem[++cnt].s=t;mem[cnt].t=s;mem[cnt].v=0;mem[cnt].next=head[t];head[t]=mem+cnt;
}
int n,m;

int S,T;
int numbs[maxn];
int d[maxn];					
edge* cur[maxn],*revpath[maxn];	

void bfs()
{
	queue<int> q;
	while(!q.empty()) q.pop();
	for (int i=1;i<=n;i++) d[i]=maxn-1;				//由初始下标决定01
	d[T]=0;q.push(T);
	while(!q.empty())
	{
		int u=q.front();
		q.pop();
		for (edge* it=head[u];it;it=it->next)
		{
			edge *now=rever(it);
			if (now->v==0||d[now->s]<n) continue;
			d[now->s]=d[u]+1;
			q.push(now->s);
		}
	}
	memset(numbs,0,sizeof(numbs));
	for (int i=1;i<=n;i++) numbs[d[i]]++;			//由初始下标决定01			
}

int isap()
{
	int flow=0;
	for (int i=1;i<=n;i++) cur[i]=head[i];			//由初始下标决定01
	int u=S;
	while(d[S]<n)
	{
		if (u==T) 
		{
			int augflow=2147483647;
			for (int i=S;i!=T;i=cur[i]->t)
				augflow=min(augflow,cur[i]->v);
			for (int i=S;i!=T;i=cur[i]->t)
			{
				cur[i]->v-=augflow;
				rever(cur[i])->v+=augflow;
			}
			flow+=augflow;u=S;
		}
		edge *e;
		for (e=cur[u];e;e=e->next)
			if (e->v&&d[u]==(d[e->t]+1)) break;
		if (e)
		{
			cur[u]=e;
			revpath[e->t]=rever(e);
			u=e->t;
		}
		else
		{
			numbs[d[u]]--;
			if (numbs[d[u]]==0) break;
			cur[u]=head[u];
			int mindist=n;
			for (edge* it=head[u];it;it=it->next)
				if (it->v) mindist=min(mindist,d[it->t]);
			d[u]=mindist+1;
			numbs[d[u]]++;
			if (u!=S) u=revpath[u]->t;
		}
	}
	return flow;
}
int main()
{
	while(~scanf("%d%d",&m,&n))
	{
		cnt=-1;
		memset(head,0,sizeof(head));
		for (int i=1;i<=m;i++)
		{
			int s,e,c;
			scanf("%d%d%d",&s,&e,&c);
			add_edge(s,e,c);
		}
		S=1;T=n;n=n;
		bfs();
		printf("%d\n",isap());
	}	
	return 0;
}
