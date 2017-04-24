#include<cstdio>
#include<cstring>
using namespace std;
#define MAXN 1000001
#define LETTER 26
struct State{
	int len, link;
	int next[LETTER];
}st[MAXN * 2];
int tree[MAXN];
int cnt, last;
void init()
{
	last = cnt = 0;
	st[cnt].len = 0; st[cnt].link = -1;
	memset(st[0].next, 0, sizeof(st[0].next));
}
inline int convert(char ch){ return ch - 'a'; }
void add(char c)
{
	c = convert(c);
	int cur = ++cnt, i;
	st[cur].len = st[last].len + 1;
	memset(st[cur].next, 0, sizeof(st[cur].next));
	for (i = last; i != -1 && !st[i].next[c]; i = st[i].link)
		st[i].next[c] = cur;
	if (i == -1)st[cur].link = 0;
	else{
		int j = st[i].next[c];
		if (st[i].len + 1 == st[j].len)
			st[cur].link = j;
		else{
			int copy = ++cnt;
			st[copy].len = st[i].len + 1;
			memcpy(st[copy].next, st[j].next, sizeof(st[j].next));
			st[copy].link = st[j].link;
			for (; i != -1 && st[i].next[c] == j; i = st[i].link)
				st[i].next[c] = copy;
			st[j].link = st[cur].link = copy;
		}
	}
	last = cur;
}
void suffixAutomation(char *s)
{
	init();
	for (int i = 0; s[i]; i++)
		add(s[i]);
}
void suffixTree(char *s)
{
	init();
	for (int i = strlen(s) - 1; i >= 0; i--){
		add(s[i]);
		tree[i] = last;
	}
}