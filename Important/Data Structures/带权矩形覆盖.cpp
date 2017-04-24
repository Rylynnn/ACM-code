//每组数据第一行n,w,h，表示星星个数(不超过60000)、天窗长度和宽度（1到1000000之间）。
//接下来n行，每行3个数，分别为x,y坐标（0<=x,y<2^31）和亮度（1到100）。
//注意，窗户不可旋转。
//输出：能看到的最大总亮度。
int n, w, h, yNum;
struct Point{
	int x, y, light, index;
	bool operator < (const Point& p)const{
		return x < p.x;
	}
}p[60001];
int y[60001], reachL[60001];
int main()
{
	while (scanf("%d%d%d", &n, &w, &h) == 3){
		int ans = 0;
		for (int i = 0; i < n; i++){
			scanf("%d%d%d", &p[i].x, &p[i].y, &p[i].light);
			y[i] = p[i].y;
		}
		sort(y, y + n);
		yNum = unique(y, y + n) - y;
		for (int i = 0, j = 0; i < yNum; i++){
			while (y[i] - y[j] >= h)j++;
			reachL[i] = j;
		}
		for (int i = 0; i < n; i++)
			p[i].index = lower_bound(y, y + yNum, p[i].y) - y;
		sort(p, p + n);
		init(n);
		for (int i = 0, j = 0; i < n; i++){
			for (; p[i].x - p[j].x >= w; j++)
				addValue(reachL[p[j].index], p[j].index + 1, p[j].light);
			addValue(reachL[p[i].index], p[i].index + 1, -p[i].light);
			ans = min(ans, tree[1].value);
		}
		printf("%d\n", -ans);
	}
}

