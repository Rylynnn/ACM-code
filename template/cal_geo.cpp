#include<bits/stdc++.h>

using namespace std;

const double eps=1e-10;
const double pi=3.1415926535897932384626433832795;
const double eln=2.718281828459045235360287471352;
//const long long RAND_MAX=1e20;
//-----------------------------------------------------------------------------------------------------------------
int dcmp(double x)
{
	if(fabs(x)<eps)return 0;else return x<0?-1:1;
}

struct point
{
	double x,y;
	point(double x=0,double y=0):x(x),y(y){}
	bool operator < (const point& b)const
	{
		return x<b.x||(x==b.x&&y<b.y);
	}
	bool operator == (const point& b)const
	{
		return dcmp(x-b.x)==0&&dcmp(y-b.y)==0;
	}
};
struct point operator + (point a,point b){return point(a.x+b.x,a.y+b.y);}
struct point operator - (point a,point b){return point(a.x-b.x,a.y-b.y);}
struct point operator * (point a,double p){return point(a.x*p,a.y*p);}
struct point operator / (point a,double p){return point(a.x/p,a.y/p);}

double Dot(point a,point b){return a.x*b.x+a.y*b.y;}
double Length(point a){return sqrt(Dot(a,a));}
double Angle(point a,point b){return acos(Dot(a,b)/(Length(a)*Length(b)));}//两向量夹角
double angle(point a){return atan2(a.y,a.x);}//向量极角
double Cross(point a,point b){return a.x*b.y-a.y*b.x;}
double Area2(point a,point b,point c){return Cross(b-a,c-a);}//三角形有向面积的两倍
double ToRad(double x){return x/180.0*pi;}
double ToDegree(double x){return x/pi*180.0;}

struct point Normal(point p)
{
	if(dcmp(Length(p))==0)return p;else return p/Length(p);
}

struct point Rotate(point a,double rad)
{
	return point(a.x*cos(rad)-a.y*sin(rad),a.x*sin(rad)+a.y*cos(rad));
}

struct point _GetLinePoint(point p,point p1,point q,point q1)//点+方向确定
{
	point u=p-q;
	double t=Cross(q1,u)/Cross(u,p1);
	return (p+p1*t);
}

struct point GetLinePoint(point p,point p1,point q,point q1)//两点确定
{
	return _GetLinePoint(p,p1-p,q,q1-q);
}

double DistanceToLine(point p,point a,point b)
{
	point v1=b-a,v2=p-a;
	return fabs(Cross(v1,v2))/Length(v1);
}

double DistanceToSegment(point p,point a,point b)
{
	if(a==b)return Length(p-a);
	point v1=b-a,v2=p-a,v3=p-b;
	if(dcmp(Dot(v1,v2))<0)return Length(v2);
	else if(dcmp(Dot(v1,v3))>0)return Length(v3);
	else return fabs(Cross(v1,v2))/Length(v1);
}

struct point PointShadow(point p,point a,point b)
{
	point v=b-a;
	return a+v*(Dot(v,p-a)/Dot(v,v));
}

bool SegmentCross(point a1,point a2,point b1,point b2)
{
	double c1=Cross(a2-a1,b1-a1);
	double c2=Cross(a2-a1,b2-a1);
	double c3=Cross(b2-b1,a1-b1);
	double c4=Cross(b2-b1,a2-b1);
	return (dcmp(c1)*dcmp(c2)<0)&&(dcmp(c3)*dcmp(c4)<0);
}

bool OS1(point p,point a,point b)
{
	double max=a.x>b.x?a.x:b.x;
	double min=a.x<b.x?a.x:b.x;
	double max1=a.y>b.y?a.y:b.y;
	double min1=a.y<b.y?a.y:b.y;
	if(p.x>=min&&p.x<=max&&p.y>=min1&&p.y<=max1)return true;else return false;
}

bool OnSegment(point p,point a,point b)
{
	if(OS1(p,a,b)&&(dcmp(Cross(a-p,b-p))==0))return true;else return false;
}

double ManySideArea(point *p,int n)//0-n-1
{
	double area=0;
	for(int i=1;i<n-1;i++)
		area+=fabs(Cross(p[i]-p[0],p[i+1]-p[0]));
	return area/2.0;
}
//--------------------------------------圆相关---------------------------------------------
struct line//p,q为确定线的两点，v为方向向量(可为单位向量也可不)
{
	point p,q,v;
	//line(point p,point v):p(p),v(v){}
	point CalPoint(double t){return p+v*t;}
};

struct circle
{
	struct point c;
	double r;
	//circle(point c,double r):c(c),r(r){}
	point Point(double a)
	{
		return point(c.x+cos(a)*r,c.y+sin(a)*r);
	}
};

int CircleAndLine(line l,circle cc,vector<point> &sol)//sol存储交点
{
	l.v=l.q-l.p;
	double a=l.v.x,b=l.p.x-cc.c.x,c=l.v.y,d=l.p.y-cc.c.y;
	double e=a*a+c*c,f=2*(a*b+c*d),g=b*b+d*d-cc.r*cc.r;
	double delta=f*f-4*e*g;
	if(dcmp(delta)<0)return 0;
	double t1,t2;
	if(dcmp(delta)==0)
	{
		t1=t2=-f/(2.0*e);sol.push_back(l.CalPoint(t1));
		return 1;
	}
	t1=(-f-sqrt(delta))/(2.0*e);sol.push_back(l.CalPoint(t1));
	t2=(-f+sqrt(delta))/(2.0*e);sol.push_back(l.CalPoint(t2));
	return 2;
}

int CircleAndCircle(circle c1,circle c2,vector<point> &sol)//sol存储交点
{
	double d=Length(c1.c-c2.c);
	if(dcmp(d)==0)
	{
		if(dcmp(c1.r-c2.r)==0)return -1;else return 0;
	}
	if(dcmp(c1.r+c2.r-d)<0)return 0;
	if(dcmp(fabs(c1.r-c2.r)-d)>0)return 0;
	
	double a=angle(c2.c-c1.c);
	double da=acos((c1.r*c1.r+d*d-c2.r*c2.r)/(2.0*c1.r*d));
	point p1=c1.Point(a-da),p2=c1.Point(a+da);

	sol.push_back(p1);
	if(p1==p2)return 1;
	sol.push_back(p2);
	return 2;
}

int GetTangent(point p,circle c,vector<line> &sol)//求过某点作圆的切线，sol存储线的信息
{//(函数名没有s)
	point u=c.c-p;
	double dist=Length(u);
	if(dcmp(dist-c.r)<0)return 0;
	else if(dcmp(dist-c.r)==0)
	{
		line l1;
		l1.p=p;
		l1.v=Normal(Rotate(u,pi/2.0));
		l1.q=l1.p+l1.v;
		sol.push_back(l1);
		return 1;
	}else
	{
		double ang=asin(c.r/dist);
		double tmp=Length(u);
		tmp=sqrt(tmp*tmp-c.r*c.r);
		line l1,l2;
		l1.p=l2.p=p;
		l1.v=Normal(Rotate(u,-ang));
		l1.q=l1.p+l1.v*tmp;sol.push_back(l1);
		l2.v=Normal(Rotate(u,+ang));
		l2.q=l2.p+l2.v*tmp;sol.push_back(l2);
		return 2;
	}
}

void GetTangentscl(bool z,line &l)
{
	if(!z)return;
	point tmp;
	tmp=l.p;l.p=l.q;l.q=tmp;
}

int GetTangents(circle a,circle b,vector<line> &sol)//两圆公切线，sol存储切线信息(p存储圆a上的切点，q存储圆b上的切点)
{//(函数名有s)
	int cnt=0;
	bool z=false;
	if(dcmp(a.r-b.r)<0)
	{
		z=true;
		circle tmp=a;a=b;b=tmp;
	}
	double d2=(a.c.x-b.c.x)*(a.c.x-b.c.x)+(a.c.y-b.c.y)*(a.c.y-b.c.y);
	double rdiff=fabs(a.r-b.r);
	double rsum=a.r+b.r;
	if(dcmp(d2-rdiff*rdiff)<0)return 0;
	
	double base=angle(b.c-a.c);
	if(dcmp(d2)==0&&dcmp(a.r-b.r)==0)return -1;//无限多条切线
	if(dcmp(d2-rdiff*rdiff)==0)//两圆内切
	{
		line l;
		point jd=a.Point(base);
		l.p=l.q=jd;
		l.v=Normal(Rotate(a.c-jd,pi/2.0));
		sol.push_back(l);return 1;
	}
	//有外公切线
	double ang=acos((a.r-b.r)/sqrt(d2));
	line l;
	l.p=a.Point(base+ang);l.q=b.Point(base+ang);GetTangentscl(z,l);l.v=Normal(l.q-l.p);sol.push_back(l);cnt++;
	l.p=a.Point(base-ang);l.q=b.Point(base-ang);GetTangentscl(z,l);l.v=Normal(l.q-l.p);sol.push_back(l);cnt++;
	if(dcmp(d2-rsum*rsum)==0)//一条内公切线
	{
		l.p=a.Point(base);l.q=b.Point(pi+base);GetTangentscl(z,l);l.v=Normal(l.q-l.p);sol.push_back(l);cnt++;
	}else if(dcmp(d2-rsum*rsum)>0)//两条内公切线
	{
		ang=acos((a.r+b.r)/sqrt(d2));
		l.p=a.Point(base+ang);l.q=b.Point(pi+base+ang);GetTangentscl(z,l);l.v=Normal(l.q-l.p);sol.push_back(l);cnt++;
		l.p=a.Point(base-ang);l.q=b.Point(pi+base-ang);GetTangentscl(z,l);l.v=Normal(l.q-l.p);sol.push_back(l);cnt++;
	}
	return cnt;
}

#define polygon vector<point> //点数组作为多边形定义

int IsPointInPolygon(point p,polygon poly)
{
	int wn=0;
	int n=poly.size();
	for(int i=0;i<n;i++)
	{
		if(OnSegment(p,poly[i],poly[(i+1)%n]))return -1;
		int k=dcmp(Cross(poly[(i+1)%n]-poly[i],p-poly[i]));
		int d1=dcmp(poly[i].y-p.y);
		int d2=dcmp(poly[(i+1)%n].y-p.y);
		if(k>0&&d1<=0&&d2>0)wn++;
		if(k<0&&d2<=0&&d1>0)wn--;
	}
	if(wn!=0)return 1;
	return 0;
}

struct dirline//有向直线类
{
	point p;
	point v;
	double ang;
	dirline(){}
	dirline(point p,point v):p(p),v(v){ang=atan2(v.y,v.x);}
	bool operator < (const dirline &l)const
	{return ang<l.ang;}
};

polygon CutPolygon(polygon poly,point a,point b)//有向直线a->b切割多边形，返回"左侧"，如果退化，可能会返回单点或线段(0-n-1保存点)(顺时针或逆时针)
{
	polygon newpoly;
	int n=poly.size();
	newpoly.clear();
	for(int i=0;i<n;i++)
	{
		point c=poly[i];
		point d=poly[(i+1)%n];
		if(dcmp(Cross(b-a,c-a))>=0)newpoly.push_back(c);
		if(dcmp(Cross(b-a,c-d))!=0)
		{
			point ip=GetLinePoint(a,b,c,d);
			if(OnSegment(ip,c,d))newpoly.push_back(ip);
		}
	}
	return newpoly;
}
//--------------------------------------半平面交-------------------------------------
bool OnLeft(dirline l,point p)
{
	return dcmp(Cross(l.v,p-l.p))>0;
}

point GetIntersection(dirline a,dirline b)//两有向直线交点，假设交点唯一存在
{
	point u=a.p-b.p;
	double t=Cross(b.v,u)/Cross(a.v,b.v);
	return a.p+a.v*t;
}

int HalfplaneIntersection(dirline* l,int n,point* poly)//半平面交主过程
{
	sort(l,l+n);//<运算符已定义，按极角排序

	int first,last;//双端队列
	point *p=new point[n];//p[i]=q[i]交q[i+1]
	dirline *q=new dirline[n];//双端队列
	q[first=last=0]=l[0];//初始化为只有一个半平面l[0]

	for(int i=1;i<n;i++)
	{
		while(first<last && !OnLeft(l[i],p[last-1]))last--;
		while(first<last && !OnLeft(l[i],p[first]))first++;
		q[++last]=l[i];
		if(dcmp(fabs(Cross(q[last].v,q[last-1].v)))==0)//两向量平行且同侧，取内侧一个
		{
			last--;
			if(OnLeft(q[last],l[i].p))q[last]=l[i];
		}
		if(first<last)p[last-1]=GetIntersection(q[last-1],q[last]);
	}

	while(first<last && !OnLeft(q[first],p[last-1]))last--;//剔除无用平面
	if(last-first<=1)return 0;
	p[last]=GetIntersection(q[last],q[first]);//计算首尾两半平面交点
	//复制到输出
	int m=0;
	for(int i=first;i<=last;i++)poly[m++]=p[i];
	return m;
}
//-----------------------------------三维计算几何----------------------------------------------
struct point3
{
	double x,y,z;
	point3(double x=0,double y=0,double z=0):x(x),y(y),z(z){}
	bool operator <(const point3 b)const
	{
		if(dcmp(x-b.x)<0)return true;
		if(dcmp(x-b.x)>0)return false;
		if(dcmp(y-b.y)<0)return true;
		if(dcmp(y-b.y)>0)return false;
		if(dcmp(z-b.z)<0)return true;
		return false;
	}
	bool operator ==(const point3 b)const
	{
		return (dcmp(x-b.x)==0)&&(dcmp(y-b.y)==0)&&(dcmp(z-b.z)==0);
	}
};

struct plane
{
	point3 p,n;
};

point3 operator +(point3 a,point3 b){return point3(a.x+b.x,a.y+b.y,a.z+b.z);}
point3 operator -(point3 a,point3 b){return point3(a.x-b.x,a.y-b.y,a.z-b.z);}
point3 operator *(point3 a,double p){return point3(a.x*p,a.y*p,a.z*p);}
point3 operator /(point3 a,double p){return point3(a.x/p,a.y/p,a.z/p);}

double Dot3(point3 a,point3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
double Length3(point3 a){return sqrt(Dot3(a,a));}
double Angle3(point3 a,point3 b){return acos(Dot3(a,b)/Length3(a)/Length3(b));}
double DistanceToPlane(point3 p,point3 p0,point3 n){return fabs(Dot3(p-p0,n));}//p0为平面上一点，n为法向量
point3 GetPlaneProjection(point3 p,point3 p0,point3 n){return p-n*Dot3(p-p0,n);}//点在平面上的投影
point3 Cross3(point3 a,point3 b){return point3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);}
point3 Normal3(point3 a){return a/Length3(a);}
double Area3(point3 a,point3 b,point3 c){return Length3(Cross3(b-a,c-a));}
point3 Law3(point3 a,point3 b,point3 c){return Cross3(b-a,c-a);}

plane ToLaw3(point3 a,point3 b,point3 c)
{
	plane pl;
	pl.p=a;
	pl.n=Law3(a,b,c);
	return pl;
}

point3 LineAndPlane(point3 p1,point3 p2,point3 p0,point3 n)//p1-p2直线与平面的交点(注意判断平行与重合，本函数不做判断(Dot(p2-p1,n)==0???))
{
	point3 v=p2-p1;
	double t=Dot3(n,p0-p1)/Dot3(n,p2-p1);//(假设不平行也不重合)
	return p1+v*t;
}

bool OnPlane(point3 p,plane pl)
{
	point3 u=p-pl.p;
	return dcmp(Dot3(u,pl.n))==0;
}

bool PointInTri(point3 p,point3 p0,point3 p1,point3 p2)//点是否在三角形内
{
	if(!OnPlane(p,ToLaw3(p0,p1,p2)))return false;
	double a1=Area3(p,p0,p1);
	double a2=Area3(p,p1,p2);
	double a3=Area3(p,p2,p0);
	return dcmp(a1+a2+a3-Area3(p0,p1,p2))==0;
}

bool TriSegIntersection(point3 p0,point3 p1,point3 p2,point3 a,point3 b,point3 &p)//线段与三角形是否相交，如果是，交点放在p返回
{
	plane pl=ToLaw3(p0,p1,p2);
	if(dcmp(Dot3(pl.n,b-a))==0)return false;
	else
	{
		double t=Dot3(pl.n,p0-a)/Dot3(pl.n,b-a);
		if(dcmp(t)<0||dcmp(t-1)>0)return false;
		p=a+(b-a)*t;
		return PointInTri(p,p0,p1,p2);
	}
}

double DistanceToLine(point3 p,point3 a,point3 b)
{
	point3 v1=b-a,v2=p-a;
	return Length3(Cross3(v1,v2))/Length3(v1);
}

double DistanceToSegment(point3 p,point3 a,point3 b)
{
	if(a==b)return Length3(p-a);
	point3 v1=b-a,v2=p-a,v3=p-b;
	if(dcmp(Dot3(v1,v2))<0)return Length3(v2);
	else if(dcmp(Dot3(v1,v3))>0)return Length3(v3);
	else return Length3(Cross3(v1,v2))/Length3(v1);
}

double Volume6(point3 a,point3 b,point3 c,point3 d)//ab，ac，ad的混合积，四面体abcd体积的6倍
{
	return Dot3(d-a,Cross3(b-a,c-a));
}
//---------------------------------------------三维凸包----------------------------------------
struct face
{
	int v[3];
	point3 normal(point3 *p)const{return Cross3(p[v[1]]-p[v[0]],p[v[2]]-p[v[0]]);}
	int cansee(point3* p,int i)const{return Dot3(p[i]-p[v[0]],normal(p))>0?1:0;}
};

double rand01(){return rand()/(double)RAND_MAX;}
double randeps(){return (rand01()-0.5)*eps;}
point3 add_noise(point3 p){return point3(p.x+randeps(),p.y+randeps(),p.z+randeps());}

vector<face> CH3D(point3* p,int n)//增量法求三维凸包，未考虑四点共面，须施加扰动
{
	int vis[n+5][n+5];
	memset(vis,0,sizeof(vis));
	vector<face> cur;
	cur.clear();
	cur.push_back((face){{0,1,2}});
	cur.push_back((face){{2,1,0}});
	for(int i=3;i<n;i++)
	{
		vector<face> next;
		next.clear();//计算每条边左侧的可见性
		for(int j=0,t=cur.size();j<t;j++,t=cur.size())
		{
			face& f=cur[j];
			int res=f.cansee(p,i);
			if(!res)next.push_back(f);
			for(int k=0;k<3;k++)vis[f.v[k]][f.v[(k+1)%3]]=res;
		}
		for(int j=0,t=cur.size();j<t;j++,t=cur.size())
			for(int k=0;k<3;k++)
			{
				int a=cur[j].v[k],b=cur[j].v[(k+1)%3];
				if(vis[a][b]!=vis[b][a]&&vis[a][b])next.push_back((face){{a,b,i}});
			}
		cur=next;
	}
	return cur;
}

int main()
{
	///*
	srand((unsigned)time(NULL));//三维凸包必须
	int i,k;
	circle c1,c2;
	point p,q;
	line l;
	scanf("%lf %lf",&c1.c.x,&c1.c.y);
	scanf("%lf",&c1.r);
	scanf("%lf %lf",&c2.c.x,&c2.c.y);
	scanf("%lf",&c2.r);
	l.p=p;l.q=q;
	vector<line> t;
	printf("%d\n",k=GetTangents(c1,c2,t));
	for(i=0;i<k;i++)printf("i:%d\na:%0.5lf %0.5lf\nb:%0.5lf %0.5lf\n",i+1,t[i].p.x,t[i].p.y,t[i].q.x,t[i].q.y);
	//*/
	return 0;
}
