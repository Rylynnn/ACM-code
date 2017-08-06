#include<cstdio>
#include<cmath>
#include<algorithm>
using namespace std;
constexpr double pi=acos(-1.0);
struct cp {
	double a,b;
	cp operator +(const cp &o) const noexcept {return (cp){a+o.a,b+o.b};}
	cp operator -(const cp &o) const noexcept {return (cp){a-o.a,b-o.b};}
	cp operator *(const cp &o) const noexcept {return (cp){a*o.a-b*o.b,b*o.a+a*o.b};}
	cp operator *(const double &o) const noexcept {return (cp){a*o,b*o};}
	cp operator !() const noexcept {return (cp){a,-b};}
}x[131183],y[131183],z[131183],w[131183];
void fft(cp * const x, const unsigned k, const unsigned v) noexcept {
	for(unsigned i=0,j=0;i<k;i++){
		if(i>j)swap(x[i],x[j]);
		for(unsigned l=k>>1;(j^=l)<l;l>>=1);
	}
	w[0]=(cp){1,0};
	for(unsigned i=2;i<=k;i<<=1){
		cp g=(cp){cos(2*pi/i),(v?-1:1)*sin(2*pi/i)};
		for(int j=(i>>1);j>=0;j-=2)w[j]=w[j>>1];
		for(unsigned j=1;j<i>>1;j+=2)w[j]=w[j-1]*g;
		for(unsigned j=0;j<k;j+=i){
			cp *a=x+j,*b=a+(i>>1);
			for(unsigned l=0;l<i>>1;l++){
				cp o=b[l]*w[l];
				b[l]=a[l]-o;
				a[l]=a[l]+o;
			}
		}
	}
	if(v)for(unsigned i=0;i<k;i++)x[i]=(cp){x[i].a/k,x[i].b/k};
}
struct buf{
	char a[33554543],*s;
	char b[33554543],*t;
	buf() noexcept:s(a),t(b){a[fread(a,1,sizeof a,stdin)]=0;}
	~buf() noexcept {fwrite(b,1,t-b,stdout);}
	operator unsigned() noexcept {
		while (*s < 48) ++s;
		return *s++ - 48;
	}
	unsigned Input() noexcept {
		unsigned x=0;
		while(*s<48)++s;
		while(*s>32)
			x=x*10+*s++-48;
		return x;
	}
	void out(int x) noexcept {
		static char c[12];
		char*i=c;
		if(!x)*t++=48;
		else{
			while(x){
				int y=x/10;
				*i++=x-y*10+48,x=y;
			}
			while(i!=c)*t++=*--i;
		}
		*t++=10;
	}
}it;
unsigned n,m,l,K;

int main() noexcept {
	n=it.Input(),m=it.Input();
	unsigned s = 0;
	for(unsigned i=0;i<=n;i++) s += (i&1?x[i>>1].b:x[i>>1].a)=it;
	for(unsigned i=0;i<=m;i++) s += (i&1?y[i>>1].b:y[i>>1].a)=it;
	
	if (!s) {
		for(unsigned i = 0; i <= n + m; ++i) {
			it.out(0);
		}
		return 0;
	}
	
	for(K=1;K<=n+m>>1;K<<=1);
	fft(x,K,0);fft(y,K,0);
    for(unsigned i = 0;i < K / 2; ++i){
        const unsigned j=K-1&K-i;
        z[i] =
            x[i]*y[i] - 
            (x[i]-!x[j]) * (y[i]-!y[j]) * (
                    w[i] + (cp){1,0}
            ) * 0.25;
    }
    for(unsigned i = K / 2;i < K; ++i){
        const unsigned j=K-1&K-i;
        z[i] =
            x[i]*y[i] - 
            (x[i]-!x[j]) * (y[i]-!y[j]) * (
                    (cp){1,0}-w[i^K>>1]
            ) * 0.25;
    }
	fft(z,K,1);
	unsigned i = 0;
    for (unsigned i = 0; i * 2 + 1 <= n + m; ++i) {
        it.out(int(z[i].a + 0.1));
        it.out(int(z[i].b + 0.1));
    }
    if ((n + m) % 2 == 0) it.out((int)(z[n + m >> 1].a + 0.1));
	return 0;
}
