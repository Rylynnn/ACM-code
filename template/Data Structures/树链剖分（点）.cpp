#include<cstdio>
#include<cstring>
#include<vector>
using namespace std;
vector<int> v[100001];
int n, cnt, color;
int father[100001], depth[100001], top[100001], id[100001], son[100001];
struct Tree{
	int maxValue, maxId, delta;
	bool set;
}tree[1 << 18];
int treeLen;
int dfs1(int i, int fa)
{
	father[i] = fa;
	depth[i] = depth[fa] + 1;
	son[i] = 0;
	int ret = 0, maxSize = 0;
	for (unsigned int j = 0; j < v[i].size(); j++){
		int t = v[i][j];
		if (t == fa)continue;
		int size = dfs1(t, i);
		ret += size;
		if (size > maxSize){
			maxSize = size;
			son[i] = t;
		}
	}
	return ret + 1;
}
void dfs2(int i, int tp)
{
	top[i] = tp;
	id[i] = cnt++;
	if (son[i])dfs2(son[i], tp);
	for (unsigned int j = 0; j < v[i].size(); j++){
		int t = v[i][j];
		if (t != father[i] && t != son[i])dfs2(t, t);
	}
}
void init()
{
	cnt = 0; depth[1] = 0;
	dfs1(1, 1);
	dfs2(1, 1);
	for (treeLen = 1; treeLen < n; treeLen *= 2);
	memset(tree, 0, sizeof(Tree) * 2 * treeLen);
}
void pushDown(int i)
{
	if (tree[i].set){
		tree[2 * i].set = tree[2 * i + 1].set = true;
		tree[2 * i].delta = tree[2 * i + 1].delta = tree[i].delta;
		tree[i].delta = 0; tree[i].set = false;
	}
	else if (tree[i].delta){
		tree[2 * i].delta += tree[i].delta;
		tree[2 * i + 1].delta += tree[i].delta;
		tree[i].delta = 0;
	}
}
int queryL, queryR;
void addInternal(int i, int l, int len)
{
	if (queryL <= l && queryR >= l + len){
		tree[i].delta++;
		return;
	}
	len >>= 1; pushDown(i);
	int mid = l + len;
	if (mid > queryL)addInternal(2 * i, l, len);
	if (mid < queryR)addInternal(2 * i + 1, mid, len);
}
inline void addValue(int l, int r){
	queryL = l; queryR = r;
	addInternal(1, 0, treeLen);
}
int addTree(int s, int t)
{
	int ret = 0;
	int top1 = top[s], top2 = top[t];
	while (top1 != top2){
		if (depth[top1] < depth[top2]){
			addValue(id[top2], id[t] + 1);
			t = father[top2]; top2 = top[t];
		}
		else{
			addValue(id[top1], id[s] + 1);
			s = father[top1]; top1 = top[s];
		}
	}
	if (depth[s] > depth[t])swap(s, t);
	addValue(id[s], id[t] + 1);
	return ret;
}
void process(int i)
{
	if (tree[i].set){
		if (tree[i].delta > tree[i].maxValue){
			tree[i].maxValue = tree[i].delta;
			tree[i].maxId = color;
		}
		return;
	}
	pushDown(i);
	process(2 * i);
	process(2 * i + 1);
}
inline void pushDownColor(int i, int j){
	if (tree[j].maxValue < tree[i].maxValue
		|| (tree[j].maxValue == tree[i].maxValue && tree[i].maxId < tree[j].maxId)){
		tree[j].maxValue = tree[i].maxValue;
		tree[j].maxId = tree[i].maxId;
	}
}
void getAns(int i)
{
	if (i < treeLen){
		pushDownColor(i, 2 * i);
		pushDownColor(i, 2 * i + 1);
		getAns(2 * i);
		getAns(2 * i + 1);
	}
}
vector<pair<int, int>> z[100001];
int main()
{
	int m;
	while (scanf("%d%d", &n, &m) == 2 && n){
		for (int i = 1; i <= n; i++)
			v[i].clear();
		for (int i = 1; i < n; i++){
			int s, t;
			scanf("%d%d", &s, &t);
			v[s].push_back(t);
			v[t].push_back(s);
		}
		for (int i = 0; i < m; i++){
			int s, t, w;
			scanf("%d%d%d", &s, &t, &w);
			z[w].push_back(make_pair(s, t));
		}
		init();
		for (color = 1; color <= 100000; color++){
			tree[1].set = true; tree[1].delta = 0;
			for (unsigned int j = 0; j < z[color].size(); j++)
				addTree(z[color][j].first, z[color][j].second);
			process(1);
			z[color].clear();
		}
		getAns(1);
		for (int i = 1; i <= n; i++)
			printf("%d\n", tree[treeLen + id[i]].maxId);
	}
	return 0;
}