#include<cstdio>
#include<cstring>
#include<queue>
using namespace std;
#define LETTER 26
struct Trie{
	int num, fail, match;
	int next[LETTER];
}trie[500001];
int cnt;
void init(){
	cnt = 1;
	memset(trie, 0, 2 * sizeof(Trie));
	trie[0].fail = 1;
}
inline int convert(char ch){ return ch - 'a'; }
void insert(char *s)
{
	int cur = 0;
	for (int i = 0; s[i]; i++){
		int &pos = trie[cur].next[convert(s[i])];
		if (!pos){
			pos = ++cnt;
			memset(&trie[cnt], 0, sizeof(Trie));
		}
		cur = pos;
	}
	trie[cur].num++;
}
void makeFail()
{
	queue<int> q; q.push(0);
	while (!q.empty()){
		int t = q.front(); q.pop();
		for (int i = 0; i < LETTER; i++){
			int &cur = trie[t].next[i];
			if (cur){
				q.push(cur);
				trie[cur].fail = trie[trie[t].fail].next[i];
				trie[cur].match = trie[cur].num ? cur : trie[trie[cur].fail].match;
			}
			else cur = trie[trie[t].fail].next[i];
		}
	}
}
int search(char *s)
{
	int ret = 0, cur = 0;
	for (int i = 0; s[i]; i++){
		cur = trie[cur].next[convert(s[i])];
		for (int temp = trie[cur].match; temp; temp = trie[trie[temp].fail].match)
			ret += trie[temp].num;
	}
	return ret;
}