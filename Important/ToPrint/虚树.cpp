#include<cstdio>
#include<algorithm>
#include<vector>
#include<cstring>
using namespace std;
#define LL long long
int cnt;
int f[200001], g[200001];
unsigned char lg2[524289];
int st[20][400000], len;
vector<pair<int,int>> v[200001];
LL h[200001];
void dfs(int i, int father)
{
	int rank = st[0][len] = ++cnt;
	f[cnt] = i; g[i] = len++;
	for (unsigned int j = 0; j < v[i].size(); j++){
		if (v[i][j].first != father){
			h[v[i][j].first] = h[i] + v[i][j].second;
			dfs(v[i][j].first, i);
			st[0][len++] = rank;
		}
	}
}
void init(int root)
{
	len = 0; cnt = 0;
	v[0].push_back(make_pair(root, 0));
	dfs(root, root);
	for (unsigned char i = 0; (1 << i) <= len; i++)
		memset(&lg2[1 << i], i, 1 << i);
	for (int i = 1; i <= lg2[len]; i++){
		for (int k = len - (1 << i); k >= 0; k--)
			st[i][k] = min(st[i - 1][k], st[i - 1][k + (1 << (i - 1))]);
	}
}
inline int query(int i, int j){
	int pos1 = g[i], pos2 = g[j];
	if (pos1 > pos2)swap(pos1, pos2);
	int t = lg2[pos2 - pos1];
	return f[min(st[t][pos1], st[t][pos2 - (1 << t) + 1])];
}
int node[200001];
int s[200001], top;
vector<int> vt[200001];
bool used[200001];//保存虚树中哪些点是原输入点
inline bool cmp(int i, int j){ return g[i] < g[j]; }
void makeVirtualTree(int n)
{
	sort(node, node + n, cmp);
	s[top = 0] = 0; node[n] = 0;
	for (int i = 0; i <= n; i++){
		int t = query(s[top], node[i]);
		if (t != s[top]){
			while (query(s[--top], t) != s[top])
				vt[s[top]].push_back(s[top + 1]);
			vt[t].push_back(s[top + 1]);
			if (s[top] != t)s[++top] = t;
		}
		s[++top] = node[i];
		used[node[i]] = 1;
	}
}
void clearVirtualTree(int i)
{
	for (unsigned int j = 0; j < vt[i].size(); j++)
		clearVirtualTree(vt[i][j]);
	vt[i].clear(); used[i] = 0;
}
LL ans[200001];
int size[200001];
LL getTotDis(int i)
{
	LL tot = 0;
	size[i] = used[i];
	for (unsigned int j = 0; j < vt[i].size(); j++){
		tot += getTotDis(vt[i][j]) + size[vt[i][j]] * (h[vt[i][j]] - h[i]);
		size[i] += size[vt[i][j]];
	}
	return tot;
}
void solve(int i, LL res, int sz)
{
	ans[i] = res;
	for (unsigned int j = 0; j < vt[i].size(); j++)
		solve(vt[i][j], res + (h[vt[i][j]] - h[i])*(sz + size[i] - 2 * size[vt[i][j]]), sz + size[i] - size[vt[i][j]]);
}
int main()
{
	int test, n, m, x, y, z;
	scanf("%d", &test);
	while (test--){
		scanf("%d%d", &n, &m);
		for (int i = 0; i <= n; i++)v[i].clear();
		for (int i = 1; i < n; i++){
			scanf("%d%d%d", &x, &y, &z);
			v[x].push_back(make_pair(y, z));
			v[y].push_back(make_pair(x, z));
		}
		init(0);
		while (m--){
			scanf("%d", &x);
			for (int i = 0; i < x; i++)
				scanf("%d", &node[i]);
			makeVirtualTree(x);
			solve(vt[0][0], getTotDis(vt[0][0]), 0);
			LL res = 1LL << 62;
			for (int i = 0; i < x; i++)
				res = min(res, ans[node[i]]);
			printf("%lld\n", res);
			clearVirtualTree(0);
		}
	}
}