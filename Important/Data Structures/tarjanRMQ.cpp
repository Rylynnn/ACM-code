#include<cstdio>
#include<vector>
using namespace std;
int a[1000001], n;
vector<int> v[1000001];
vector<pair<int,int>> q[1000001];
int f[1000001], res[1000001];
int cartesianTree()
{
	for (int i = 0; i < n; i++)
		v[i].clear();
	vector<int> s;
	s.push_back(0);
	for (int i = 1; i < n; i++){
		int t, size = s.size();
		for (; !s.empty() && a[i]<a[s.back()]; s.pop_back())
			t = s.back();
		if (!s.empty())v[s.back()].push_back(i);
		if (size != s.size())v[i].push_back(t);
		s.push_back(i);
	}
	return s[0];
}
int getFather(int i)
{
	if (f[i] == i)return i;
	return f[i] = getFather(f[i]);
}
void tarjan(int root)
{
	f[root] = root;
	for (unsigned int i = 0; i < v[root].size(); i++){
		tarjan(v[root][i]);
		f[v[root][i]] = root;
	}
	for (unsigned int i = 0; i < q[root].size(); i++)
		res[q[root][i].second] = getFather(q[root][i].first);
}
//下标从0开始，不用清空任何数组
void query()
{
	int m;
	scanf("%d", &m);
	for (int i = 0; i < n; i++){
		q[i].clear();
		v[i].clear();
	}
	for (int i = 0; i < m; i++){
		int s, t;
		scanf("%d%d", &s, &t);
		q[t].push_back(make_pair(s,i));
	}
	tarjan(cartesianTree());
}