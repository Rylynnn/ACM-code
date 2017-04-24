#include<cstdio>
#include<cstring>
#include<algorithm>
using namespace std;
#define MAXN 50001
#define LL long long
struct Edge{
	int t, w, next;
	LL c;
}e[1000001];
int head[MAXN];
LL d[MAXN];
bool used[MAXN];
int cnt, src, dst, cur;
void init(int n)
{
	cur = 0; cnt = n;
	memset(head + 1, -1, sizeof(int)*n);
}
void addEdge(int s, int t, int w, LL c)
{
	e[cur] = { t, w, head[s], c };
	head[s] = cur++;
	e[cur] = { s, 0, head[t], -c };
	head[t] = cur++;
}
bool dijkstra()
{
	static pair<LL, int> q[MAXN];
	memset(d, 0x3f, sizeof(LL)*(cnt + 1));
	memset(used + 1, 0, sizeof(bool)*cnt);
	d[dst] = 0; q[0] = make_pair(0, dst);
	for (int pos = 1; pos;){
		int i = q->second;
		pop_heap(q, q + pos--);
		if (used[i])continue;
		used[i] = true;
		for (int j = head[i]; j != -1; j = e[j].next){
			int t = e[j].t;
			if (e[j ^ 1].w && d[t] > d[i] - e[j].c){
				d[t] = d[i] - e[j].c;
				q[pos++] = make_pair(-d[t], t);
				push_heap(q, q + pos);
			}
		}
	}
	return d[src] < d[0];
}
/*bool dijkstra()
{
	memset(d, 0x3f, sizeof(LL)*(cnt + 1));
	memset(used + 1, 0, sizeof(bool)*cnt);
	d[dst] = 0;
	for (int i = dst; d[i] != d[0];){
		used[i] = true;
		for (int j = head[i]; j != -1; j = e[j].next){
			if (e[j ^ 1].w)
				d[e[j].t] = min(d[e[j].t], d[i] - e[j].c);
		}
		i = 0;
		for (int j = 1; j <= cnt; j++){
			if (d[i] > d[j] && !used[j])i = j;
		}
	}
	return d[src] < d[0];
}*/
int dfs(int i, int flow)
{
	if (i == dst)return flow;
	used[i] = true;
	int ret = 0;
	for (int j = head[i]; j != -1; j = e[j].next){
		if (!e[j].c && e[j].w && !used[e[j].t]){
			int w = dfs(e[j].t, min(flow - ret, e[j].w));
			e[j].w -= w; e[j ^ 1].w += w; ret += w;
			if (ret == flow)break;
		}
	}
	if (ret)used[i] = false;
	return ret;
}
LL costFlow()
{
	LL ret = 0, dis = 0;
	while (dijkstra()){
		for (int i = 1; i <= cnt; i++){
			for (int j = head[i]; j != -1; j = e[j].next)
				e[j].c += d[e[j].t] - d[i];
		}
		dis += d[src];
		memset(used + 1, 0, sizeof(bool)*cnt);
		ret += dis * dfs(src, 0x7fffffff);
	}
	return ret;
}
