LL a[340000], b[340000];
int main()
{
	long long n;
	while (scanf("%lld", &n) == 1){
		int s = (int)sqrt(n + 0.5);
		for (int i = 1; i <= s; i++){
			a[i] = n / i - 1;
			b[i] = i - 1;
		}
		for (int i = 2; i <= s; i++){
			if (b[i - 1] == b[i])continue;
			LL cur = b[i - 1], t = n / i;
			int j = 1;
			for (int k = i; k <= s; j++, k += i)
				a[j] -= a[k] - cur;
			for (int k = min((LL)s, (t / i)); k >= j; k--)
				a[k] -= b[t / k] - cur;
			for (LL k = s; k >= (LL)i*i; k--)
				b[k] -= b[k / i] - cur;
		}
		printf("%lld\n", a[1]);
	}
	return 0;
}
