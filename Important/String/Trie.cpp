#include<cstdio>
#include<cstring>
#define LETTER 26
struct Trie{
	int num;
	int next[LETTER];
}trie[500001];
int cnt;
void init(){
	cnt = 0;
	memset(trie, 0, sizeof(Trie));
}
inline int convert(char ch){ return ch - 'a'; }
inline char convert2(int value){ return value + 'a'; }
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
int search(char *s)
{
	int cur = 0;
	for (int i = 0; s[i]; i++){
		cur = trie[cur].next[convert(s[i])];
		if (!cur)return -1;
	}
	return trie[cur].num;
}
char temp[1001];
void dfs(int i, int h)
{
	if (trie[i].num){
		temp[h] = 0;
		printf("%s %d\n", temp, trie[i].num);
	}
	for (int j = 0; j < LETTER; j++){
		if (trie[i].next[j]){
			temp[h] = convert2(j);
			dfs(trie[i].next[j], h + 1);
		}
	}
}