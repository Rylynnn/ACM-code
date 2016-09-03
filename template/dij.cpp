//dij heatwv
#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<algorithm>
#include<cstring>
#include<cmath>
#include<queue>
using namespace std;
#define pii pair<int, int>
priority_queue<pii, vector<pii>, greater<pii> > heap; 
struct edge
{
	int v, l, next;
}e[13007];
int n, m, S, T;
int tot=2, dist[2502], head[2502];
bool visited[2502];
void addedge(int x,int y,int z)
{
	e[tot].v=y, e[tot].l=z, e[tot].next=head[x], head[x]=tot++;
	e[tot].v=x, e[tot].l=z, e[tot].next=head[y], head[y]=tot++;
}
void dij(int x)
{
	if (x == T) return ;
	visited[x] = true;
	for (int p = head[x]; p; p = e[p].next)
		if (!visited[e[p].v] && dist[e[p].v] > (dist[x] + e[p].l))
			dist[e[p].v] = dist[x] + e[p].l, heap.push(make_pair(dist[e[p].v], e[p].v));
	while (!heap.empty() && visited[heap.top().second])
		heap.pop();
	dij(heap.top().second);
}		
int main()
{
	scanf("%d%d%d%d",&n,&m,&S,&T);
	int x, y, z;
	for (int i=1;i<=m;i++)
		scanf("%d%d%d",&x,&y,&z), addedge(x, y, z);
	for (int i=1;i<=n;i++)
		dist[i]=0x7ffffff;
	dist[S]=0;
	dij(S);
	printf("%d\n",dist[T]);
	return 0;
}
