#include<cstdio>
#include<vector>
using namespace std;
vector<pair<int, int>> v[200001];//边及该边的编号
int w[200001];//边权
int n, cnt;
int father[200001], depth[200001], top[200001], id[200001];
int f[200001];//边在树状数组（线段树）中的位置
int tmp[200001];
int dfs1(int i, int fa)
{
	father[i] = fa;
	depth[i] = depth[fa] + 1;
	tmp[i] = -1;
	int ret = 0, maxSize = 0;
	for (unsigned int j = 0; j < v[i].size(); j++){
		int t = v[i][j].first;
		if (t == fa)continue;
		int size = dfs1(t, i);
		ret += size;
		if (size > maxSize){
			maxSize = size;
			tmp[i] = j;
		}
	}
	return ret + 1;
}
void dfs2(int i, int tp, int index)
{
	top[i] = tp;
	id[i] = cnt;
	f[index] = cnt++;
	if (tmp[i] != -1)
		dfs2(v[i][tmp[i]].first, tp, v[i][tmp[i]].second);
	for (unsigned int j = 0; j < v[i].size(); j++){
		int t = v[i][j].first;
		if (t != father[i] && j != tmp[i])
			dfs2(t, t, v[i][j].second);
	}
}
int queryTree(int s, int t)
{
	int ret = 0;
	int top1 = top[s], top2 = top[t];
	while (top1 != top2){
		if (depth[top1] < depth[top2]){
			ret += sum(id[t]) - sum(id[top2] - 1);
			t = father[top2]; top2 = top[t];
		}
		else{
			ret += sum(id[s]) - sum(id[top1] - 1);
			s = father[top1]; top1 = top[s];
		}
	}
	if (s != t){
		if (depth[s] > depth[t])swap(s, t);
		ret += sum(id[t]) - sum(id[s]);
	}
	return ret;
}
void init()
{
	cnt = 0;
	dfs1(1, 1);
	dfs2(1, 1, 0);
	for (int i = 1; i < n; i++)
		tree[f[i]] = w[i];
	build();
}
int main()
{
	int q, cur;
	scanf("%d%d%d", &n, &q, &cur);
	for (int i = 1; i < n; i++){
		int s, t, value;
		scanf("%d%d%d", &s, &t, &value);
		v[s].push_back(make_pair(t, i));
		v[t].push_back(make_pair(s, i));
		w[i] = value;
	}
	init();
	for (int i = 0; i < q; i++){
		int t, u, value;
		scanf("%d%d", &t, &u);
		if (t == 0){
			printf("%d\n", queryTree(cur, u));
			cur = u;
		}
		else{
			scanf("%d", &value);
			add(f[u], value - w[u]);
			w[u] = value;
		}
	}
	return 0;
}