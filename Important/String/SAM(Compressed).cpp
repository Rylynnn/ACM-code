#include<cstdio>
#include<cstring>
using namespace std;
#define MAXN 1000001
#define LETTER 26
int pool[MAXN * LETTER / 3], poolNum;
struct State{
	int len, link;
	char ch[4];
	int next[3];
	int* child(char c){
		if (ch[0] == -1)
			return pool[next[0] + c] ? &pool[next[0] + c] : 0;
		for (int i = ch[0]; i; i--){
			if (ch[i] == c)return &next[i - 1];
		}
		return 0;
	}
	void insert(char c, int pos){
		if (ch[0] == 3){
			ch[0] = -1;
			memset(pool + poolNum, 0, sizeof(int)*LETTER);
			for (int i = 3; i; i--)
				pool[poolNum + ch[i]] = next[i - 1];
			next[0] = poolNum;
			poolNum += LETTER;
		}
		if (ch[0] == -1)pool[next[0] + c] = pos;
		else{
			next[ch[0]] = pos;
			ch[++ch[0]] = c;
		}
	}
}st[MAXN * 2];
int cnt, last;
void init()
{
	last = cnt = 0;
	st[cnt].len = 0; st[cnt].link = -1;
	st[cnt].ch[0] = 0;
	poolNum = 0;
}
inline int convert(char ch){ return ch - 'a'; }
void add(char c)
{
	c = convert(c);
	int cur = ++cnt, i, *tmp;
	st[cur].len = st[last].len + 1;
	st[cur].ch[0] = 0;
	for (i = last; i != -1 && !(tmp = st[i].child(c)); i = st[i].link)
		st[i].insert(c, cur);
	if (i == -1)st[cur].link = 0;
	else{
		int j = *tmp;
		if (st[i].len + 1 == st[j].len)
			st[cur].link = j;
		else{
			int copy = ++cnt;
			st[copy].len = st[i].len + 1;
			if (st[j].ch[0] == -1){
				st[copy].ch[0] = -1;
				st[copy].next[0] = poolNum;
				memcpy(pool + poolNum, pool + st[j].next[0], sizeof(int)*LETTER);
				poolNum += LETTER;
			}
			else{
				for (int i = 0; i < 4; i++)
					st[copy].ch[i] = st[j].ch[i];
				for (int i = st[j].ch[0]; i; i--)
					st[copy].next[i - 1] = st[j].next[i - 1];
			}
			st[copy].link = st[j].link;
			for (; i != -1 && (tmp = st[i].child(c)) && *tmp == j; i = st[i].link)
				*tmp = copy;
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