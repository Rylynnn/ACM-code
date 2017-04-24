#include<algorithm>
using namespace std;
int maxCommonSubstring(char *s1, char *s2)
{
	init();
	suffixAutomation(s1);
	int match = 0, ret = 0, cur = 0;
	for (int i = 0; s2[i]; i++){
		char c = convert(s2[i]);
		if (!st[cur].next[c]){
			while (cur != -1 && !st[cur].next[c])
				cur = st[cur].link;
			if (cur == -1){ match = cur = 0; continue; }
			match = st[cur].len;
		}
		cur = st[cur].next[c];
		ret = max(ret, ++match);
	}
	return ret;
}