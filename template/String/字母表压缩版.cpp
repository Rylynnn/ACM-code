#define LETTER 26
struct Trie{
	int num, next, fail;
}trie[1000000];
int cnt;
int pool[LETTER * 200000], poolEnd;
void init()
{
	cnt = 0;
	trie[0].num = 0;
	trie[0].next = -1;
	memset(pool, 0, sizeof(pool));
	poolEnd = 0;
}
inline int convert(char ch){ return ch - 'a'; }
inline bool oneBranch(int value){ return value < LETTER; }
inline int child(int i, int ch){
	if (oneBranch(trie[i].next))return trie[i].next == ch ? i + 1 : 0;
	return pool[trie[i].next + ch];
}
void insert(char *s)
{
	int pos = 0, i;
	for (i = 0; s[i]; i++){
		int t = trie[pos].next;
		if (oneBranch(t)){
			if (t == convert(s[i]))pos++;
			else{
				trie[pos].next = (poolEnd += LETTER);
				if (t != -1)pool[trie[pos].next + t] = pos + 1;
				break;
			}
		}
		else if (pool[t + convert(s[i])])
			pos = pool[t + convert(s[i])];
		else break;
	}
	if (s[i]){
		pool[trie[pos].next + convert(s[i])] = ++cnt;
		for (i++; s[i]; i++, cnt++){
			trie[cnt].num = 0;
			trie[cnt].next = convert(s[i]);
		}
		trie[cnt].num = 1;
		trie[cnt].next = -1;
	}
	else trie[pos].num++;
}
int getFailPoint(int father, int ch)
{
	while (father){
		father = trie[father].fail;
		int pos = child(father, ch);
		if (pos)return pos;
	}
	return 0;
}
void makeFail()
{
	queue<int> q; q.push(0);
	trie[0].fail = 0;
	while (!q.empty()){
		int t = q.front(); q.pop();
		if (oneBranch(trie[t].next)){
			if (trie[t].next != -1){
				trie[t + 1].fail = getFailPoint(t, trie[t].next);
				q.push(t + 1);
			}
		}
		else for (int i = 0; i < LETTER; i++){
			int cur = pool[trie[t].next + i];
			if (cur){
				trie[cur].fail=getFailPoint(t, i);
				q.push(cur);
			}
		}
	}
}
//统计匹配总次数，包括母串多次匹配同一模式串或多个模式串相同
int search(char *s)
{
	int ret = 0, cur = 0;
	for (int i = 0; s[i]; i++){
		int ch = convert(s[i]);
		for (; cur && !child(cur, ch); cur = trie[cur].fail);
		cur = child(cur, ch);
		for (int temp = cur; temp; temp = trie[temp].fail)
			ret += trie[temp].num;
	}
	return ret;
}
