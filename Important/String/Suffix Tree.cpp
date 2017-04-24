#include<vector>
using namespace std;
int bucket[MAXN], order[MAXN];
int id[MAXN];
bool leave[MAXN];
vector<int> tree[MAXN];
void sortState()//后缀树层序遍历
{
	int size = st[last].len;
	memset(bucket, 0, sizeof(int)*(size + 1));
	for (int i = 0; i <= cnt; i++)
		bucket[st[i].len]++;
	for (int i = 1; i <= size; i++)
		bucket[i] += bucket[i - 1];
	for (int i = cnt; i; i--)
		order[--bucket[st[i].len]] = i;
}
void suffixTree(char *s)
{
	int len = strlen(s);
	init();
	memset(leave, 0, len * sizeof(bool) * 2);
	for (int i = len - 1; i >= 0; i--){
		add(s[i]);
		id[i] = last;
		leave[last] = true;
	}
	for (int i = 0; i <= cnt; i++)
		tree[i].clear();
	for (int i = cnt; i; i--)
		tree[st[i].link].push_back(i);
}