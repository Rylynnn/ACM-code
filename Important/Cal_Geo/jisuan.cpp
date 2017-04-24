#include<iostream>
#include<cmath>
#include<algorithm>
#include<vector>
using namespace std;
#define EPS 1e-10
#define INF 1e10
#define PI 3.14159265358979323846
inline bool EQUAL(double t1, double t2){
	return t1 - t2 < EPS && t1 - t2 > -EPS;
}
inline bool LESS(double t1, double t2){
	return t1 <= t2 - EPS;
}
inline bool LESS_EQUAL(double t1, double t2){
	return t1 < t2 + EPS;
}
inline int SGN(double t){
	return LESS(t, 0) ? -1 : LESS(0, t) ? 1 : 0;
}
class Point
{
public:
	double x, y;
	Point(){}
	Point(double x, double y) :x(x), y(y){}

	bool operator == (const Point& p)const{
		return EQUAL(x, p.x) && EQUAL(y, p.y);
	}
	bool operator < (const Point& p)const{
		return LESS_EQUAL(x, p.x) && (LESS(x, p.x) || LESS(y, p.y));
	}
	Point operator + (const Point& p)const{
		return Point(x + p.x, y + p.y);
	}
	Point operator - (const Point& p)const{
		return Point(x - p.x, y - p.y);
	}
	double operator * (const Point& p)const{
		return x*p.y - y*p.x;
	}
	Point operator * (double value)const{
		return Point(x*value, y*value);
	}
	Point operator / (double value)const{
		return Point(x / value, y / value);
	}
	double dot(const Point& p)const{
		return x*p.x + y*p.y;
	}
	double r2()const{ return x*x + y*y; }
	double r()const{ return hypot(x, y); }
	double dis2(const Point& p)const{
		return (*this - p).r2();
	}
	double dis(const Point& p)const{
		return (*this - p).r();
	}

	bool onLine(const Point& p1, const Point& p2)const{
		return EQUAL((*this - p1)*(*this - p2), 0);
	}
	bool onLineSeg(const Point& p1, const Point& p2)const{
		//include extream points
		return onLine(p1, p2) && inRect(p1, p2);
	}
	double lineRelation(const Point& p1, const Point& p2)const{
		Point t = p2 - p1;
		return t.dot(*this - p1) / t.r2();
		//ret 0, *this=p1; ret 1,*this=p2;
		//ret (0,1), *this is interior to p1p2
	}
	Point footPoint(const Point& p1, const Point& p2)const{
		double r = lineRelation(p1, p2);
		return p1 + (p2 - p1)*r;
	}
	double lineDis(const Point& p1, const Point& p2)const{
		return abs((p1 - *this)*(p2 - *this)) / p1.dis(p2);
	}
	double lineSegDis(const Point& p1, const Point& p2, Point& ret)const;
	double lineSegArrayDis(const Point* p, int lineNum, Point& ret)const;
	bool lineSegArrayDisCmp(const Point* p, int lineNum, double value)const;
	Point mirror(Point& p1, Point& p2){
		Point foot = footPoint(p1, p2);
		return foot * 2 - *this;
	}

	Point rotate(double angle)const{
		Point f(sin(angle), cos(angle));
		return Point(*this * f, dot(f));
	}
	Point rotate90()const{
		return Point(-y, x);
	}
	double cosAngle(const Point& p1, const Point& p2)const{
		Point t1 = *this - p1, t2 = *this - p2;
		return t1.dot(t2) / sqrt(t1.r2()*t2.r2());
	}
	double tanAngle(const Point& o = Point(0, 0))const{
		if (EQUAL(x, o.x))return y - o.y >= 0 ? INF : -INF;
		return (y - o.y) / (x - o.x);
	}
	double angle(const Point& p1, const Point& p2)const{
		return acos(cosAngle(p1, p2));
	}
	double angle(const Point& o = Point(0, 0))const{
		return atan2(y - o.y, x - o.x);
	}
	//left return 1, right return -1, on line return 0.
	int direction(const Point& p1, const Point& p2)const{
		return SGN(x*(p1.y - p2.y) + p1.x*(p2.y - y) + p2.x*(y - p1.y));
	}
	
	bool inRect(const Point& p1, const Point& p2)const{
		return LESS_EQUAL((p1.x - x)*(p2.x - x), 0) && LESS_EQUAL((p1.y - y)*(p2.y - y), 0);
	}
	int inPolygon(const Point* p, int n)const;
	int inConvex(const Point* p, int n)const;
	int inCircle(const Point& o, double r)const{
		double dist = dis2(o);
		return SGN(r*r - dist);
	}
	void pointcut(const Point& o, double r, Point& ret1, Point& ret2)const;
	Point nearnestPoint(const Point& o, double r)const;
};
double Point::lineSegDis(const Point& p1, const Point& p2, Point& ret)const
{
	double r = lineRelation(p1, p2);
	if (LESS_EQUAL(r, 0))ret = p1;
	else if (LESS_EQUAL(1, r))ret = p2;
	else ret = footPoint(p1, p2);
	return dis(ret);
}
//input lineNum+1 points
double Point::lineSegArrayDis(const Point* p, int lineNum, Point& ret)const
{
	Point tp;
	double td, mind = INF;
	for (int i = 0; i < lineNum; i++){
		td = lineSegDis(p[i], p[i + 1], tp);
		if (LESS(td, mind)){
			mind = td; ret = tp;
		}
	}
	return mind;
}
//input lineNum+1 points
bool Point::lineSegArrayDisCmp(const Point* p, int lineNum, double value)const
{
	Point tp;
	double td;
	int flag = 1;
	for (int i = 0; i < lineNum; i++){
		td = lineSegDis(p[i], p[i + 1], tp);
		if (LESS_EQUAL(td, value))
			return true;
	}
	return false;
}

//donnot include extream points, and donnot include coincidence.
inline bool lineSegLineSegIntersect(const Point& p1, const Point& p2, const Point& q1, const Point& q2){
	Point pq1 = p1 - q1, p12 = p2 - p1, q12 = q2 - q1;
	return SGN(pq1*q12)*SGN((p2 - q1)*q12) < 0 && SGN(pq1*p12)*SGN((p1 - q2)*p12) < 0;
}
//include extream points and coincidence.
inline bool lineSegLineSegIntersect2(const Point& p1, const Point& p2, const Point& q1, const Point& q2){
	if (!(LESS_EQUAL(min(q1.x, q2.x), max(p1.x, p2.x)) && LESS_EQUAL(min(p1.x, p2.x), max(q1.x, q2.x))
		&& LESS_EQUAL(min(q1.y, q2.y), max(p1.y, p2.y)) && LESS_EQUAL(min(p1.y, p2.y), max(q1.y, q2.y))))
		return false;
	Point pq1 = p1 - q1, p12 = p2 - p1, q12 = q2 - q1;
	return SGN(pq1*q12)*SGN((p2 - q1)*q12) <= 0 && SGN(pq1*p12)*SGN((p1 - q2)*p12) <= 0;
}
//donot include extream points, and donot include coincidence.
inline bool lineLineSegIntersect(const Point& l1, const Point& l2, const Point& p1, const Point& p2){
	Point line = l2 - l1;
	return SGN((p1 - l1)*line)*SGN((p2 - l1)*line) < 0;
}
//donnot include coincidence.
inline bool lineLineIntersect(const Point& p1, const Point& p2, const Point& q1, const Point& q2){
	return !EQUAL((p2 - p1)*(q2 - q1), 0);
}
inline Point lineLineIntersectPoint(const Point& p1, const Point& p2, const Point& q1, const Point& q2){
	Point q12 = q2 - q1;
	double k = (p2 - p1)*q12;
	if (EQUAL(k, 0))return Point(INF*INF, INF*INF);
	double r = ((q1 - p1)*q12) / k;
	return p1 + (p2 - p1) * r;
}

Point circumcenter(const Point& p1, const Point& p2, const Point& p3)
{
	Point t1 = (p1 + p2)*0.5, t2, t3 = (p2 + p3)*0.5, t4;
	t2 = t1 + (p1 - p2).rotate90();
	t4 = t3 + (p2 - p3).rotate90();
	return lineLineIntersectPoint(t1, t2, t3, t4);
}
Point incenter(const Point& p1, const Point& p2, const Point& p3)
{
	double r12 = p1.dis(p2), r23 = p2.dis(p3), r31 = p3.dis(p1);
	Point t1 = (p2*r31 + p3*r12) / (r12 + r31), t2 = (p1*r23 + p3*r12) / (r12 + r23);
	return lineLineIntersectPoint(p1, t1, p2, t2);
}
Point prepencenter(const Point& p1, const Point& p2, const Point& p3)
{
	Point t1 = p1 + (p2 - p3).rotate90();
	Point t2 = p2 + (p1 - p3).rotate90();
	return lineLineIntersectPoint(p1, t1, p2, t2);
}
inline Point barycenter(const Point& p1, const Point& p2, const Point& p3){
	return (p1 + p2 + p3) / 3;
}
inline double apothem(const Point& p1, const Point& p2, const Point& p3){
	Point p12 = p2 - p1, p13 = p3 - p1, p23 = p3 - p2;
	return abs(p12*p23) / (p12.r() + p13.r() + p23.r());
}
inline double circumradius(const Point& p1, const Point& p2, const Point& p3){
	Point p12 = p2 - p1, p13 = p3 - p1, p23 = p3 - p2;
	return sqrt(p12.r2()*p23.r2()*p13.r2()) / (2 * abs(p12*p23));
}

int getPolygonDirection(const Point* p, int n)
{
	int index = 0;
	for (int i = 1; i < n; i++){
		if (p[i] < p[index])index = i;
	}
	return p[index].direction(p[index + 1 < n ? index + 1 : 0], p[index - 1 >= 0 ? index - 1 : n - 1]);
}
bool checkConvex(const Point* p, int n)
{
	int direction = p[0].direction(p[n - 1], p[1]);
	if (direction == 0)return false;
	if (p[n - 1].direction(p[n - 2], p[0]) != direction)return false;
	for (int i = n - 2; i > 0; i--){
		if (p[i].direction(p[i - 1], p[i + 1]) != direction)
			return false;
	}
	return true;
}
bool checkConvex(const Point* p, int n, bool *ret)
{
	bool retValue = true;
	int direction = getPolygonDirection(p, n);
	if (!(ret[n - 1] = p[n - 1].direction(p[0], p[n - 2]) == direction))
		retValue = false;
	if (!(ret[0] = p[0].direction(p[1], p[n - 1]) == direction))
		retValue = false;
	for (int i = n - 2; i > 0; i--){
		if (!(ret[i] = p[i].direction(p[i + 1], p[i - 1]) == direction))
			retValue = false;
	}
	return retValue;
}
double polygonArea(const Point* p, int n)
{
	double area = 0;
	for (int i = n - 2; i > 0; i--)
		area += p[i].y *(p[i - 1].x - p[i + 1].x);
	area += p[0].y*(p[n - 1].x - p[1].x);
	area += p[n - 1].y*(p[n - 2].x - p[0].x);
	return area / 2;
}
int Point::inPolygon(const Point* p, int n)const
{
	int i, j = n - 1, odd = -1;
	for (i = 0; i < n; j = i++){
		if (LESS(p[i].y, y) != LESS(p[j].y, y)){
			double tx = (y - p[j].y) / (p[i].y - p[j].y)*(p[i].x - p[j].x) + p[j].x;
			if (LESS_EQUAL(tx, x)){
				if (LESS(tx, x))odd = -odd;
				else return 0;
			}
		}
		else if (onLineSeg(p[i], p[j]))return 0;
	}
	return odd;
}
int Point::inConvex(const Point* p, int n)const
{
	int _direction = p[1].direction(p[2], p[0]);
	if (direction(p[0], p[1]) != _direction){
		if (onLineSeg(p[0], p[1]))return 0;
		return -1;
	}
	if (direction(p[n - 1], p[0]) != _direction){
		if (onLineSeg(p[n - 1], p[0]))return 0;
		return -1;
	}
	int left = 2, right = n - 1;
	while (left < right){
		int mid = (left + right) >> 1;
		if (direction(p[0], p[mid]) == _direction)left = mid + 1;
		else right = mid;
	}
	int ret = direction(p[left-1],p[left]);
	return ret == _direction ? 1 : ret == 0 ? 0 : -1;
}
Point lineConvexIntersectPointInternal(const Point& p1, const Point& p2, const Point* p, int n, int start, int end)
{
	Point p12 = p2 - p1;
	if (end < start)end += n;
	double value = SGN((p[start] - p1)*p12);
	while (start + 1 < end){
		int mid = (start + end) / 2;
		Point cur = p[mid < n ? mid : mid - n];
		double t = (cur - p1)*p12*value;
		if (LESS(0, t))start = mid;
		else if (LESS(t, 0))end = mid;
		else return cur;
	}
	if (start >= n)start -= n;
	return lineLineIntersectPoint(p1, p2, p[start], p[start + 1]);
}
int lineConvexIntersectPoint(const Point& p1, const Point& p2, const Point* p, int n, Point& ret1, Point& ret2)
{
	Point p12 = p2 - p1;
	int pos = 0, step = n * 2 / 3;
	double d = (p[pos] - p1)*p12;
	int zero = -1, pos2 = -1;
	while (step > 1){
		step=(step + 1) / 2;
		int i = pos + step, k = pos - step;
		if (i >= n)i -= n;
		if (k < 0)k += n;
		double di = (p[i] - p1)*p12, dk = (p[k] - p1)*p12;
		if (SGN(di)*SGN(d) < 0){ pos2 = i; break; }
		if (SGN(dk)*SGN(d) < 0){ pos2 = k; break; }
		if (abs(di) < abs(d)){ d = di; pos = i; }
		if (abs(dk) < abs(d)){ d = dk; pos = k; }
		if (EQUAL(d, 0)){ zero = pos; break; }
	}
	if (zero != -1){
		ret1 = p[zero];
		int left = zero - 1 >= 0 ? zero - 1 : n - 1;
		int right = zero + 1 < n ? zero + 1 : 0;
		double dl = (p[left] - p1)*p12, dr = (p[right] - p1)*p12;
		if (EQUAL(dl, 0)){ ret2 = p[left]; return 3; }
		else if (EQUAL(dr, 0)){ ret2 = p[right]; return 3; }
		else if (dl*dr < 0)return 1;
		else{ pos = left; pos2 = right; }
	}
	if (pos2 == -1)return 0;
	ret1 = lineConvexIntersectPointInternal(p1, p2, p, n, pos, pos2);
	ret2 = lineConvexIntersectPointInternal(p1, p2, p, n, pos2, pos);
	return 2;
}

bool lineSegInPolygon(const Point& p1, const Point& p2, const Point* p, int n)
{
	bool flag = false;
	Point minPoint;
	switch (p1.inPolygon(p, n)){
	case -1:return false;
	case 0:flag = true;
	}
	switch (p2.inPolygon(p, n)){
	case -1:return false;
	case 1:flag = false;
	}
	if (flag)minPoint = max(p1, p2);
	for (int i = 0, j = n - 1; i < n; j = i++){
		if (p[i].onLineSeg(p1, p2) && !(p[i] == p1 || p[i] == p2)){
			if (p[i > 0 ? i - 1 : n - 1].direction(p1, p2) * p[i + 1 < n ? i + 1 : 0].direction(p1, p2) < 0)
				return false;
			if (flag && p[i] < minPoint)minPoint = p[i];
		}
		else if (lineSegLineSegIntersect(p[i], p[j], p1, p2))
			return false;
	}
	if (flag){
		const Point& t = min(p1, p2);
		Point mid = (t + minPoint)*0.5;
		if (mid.inPolygon(p, n) == -1)return false;
	}
	return true;
}
Point gravityCenter(const Point* p, int n)
{
	if (n < 3){
		if (n == 1)return p[0];
		else return (p[0] + p[1])*0.5;
	}
	double area = 0;
	Point ret(0, 0);
	for (int i = 0, j = n - 1; i < n; j = i++){
		double t = p[i] * p[j];
		area += t;
		ret.x += (p[i].x + p[j].x)*t;
		ret.y += (p[i].y + p[j].y)*t;
	}
	return ret / (3 * area);
}
//ret[n] must be available to visit.
int convexHullSorted(const Point* p, int n, Point* ret)
{
	int j = 0;
	for (int i = 0; i < n; i++){
		while (j >= 2 && p[i].direction(ret[j - 2], ret[j - 1]) != 1)j--;
		ret[j++] = p[i];
	}
	int mid = j + 1;
	for (int i = n - 2; i >= 0; i--){
		while (j >= mid && p[i].direction(ret[j - 2], ret[j - 1]) != 1)j--;
		ret[j++] = p[i];
	}
	return j - 1;
}
void convexHullSorted(const Point* p, int n, Point* up, int& retUp, Point* down, int& retDown)
{
	retUp = retDown = 0;
	for (int i = 0; i < n; i++){
		while (retUp >= 2 && p[i].direction(up[retUp - 2], up[retUp - 1]) != -1)retUp--;
		while (retDown >= 2 && p[i].direction(down[retDown - 2], down[retDown - 1]) != 1)retDown--;
		up[retUp++] = p[i];
		down[retDown++] = p[i];
	}
}
int halfPlainIntersectInternal(vector<pair<double, const Point*>>& v, int n, Point* ret)
{
	for (int i = 0; i < n; i++)
		v[i].first = v[i].second[1].angle(v[i].second[0]);
	sort(v.begin(), v.end());
	vector<const Point*> line(n);
	vector<Point> point(n);
	int first = 0, last = 0;
	line[0] = v[0].second;
	for (unsigned int i = 1; i < v.size(); i++){
		while (first < last && point[last - 1].direction(v[i].second[0], v[i].second[1]) == -1) last--;
		while (first < last && point[first].direction(v[i].second[0], v[i].second[1]) == -1) first++;
		line[++last] = v[i].second;
		if (!lineLineIntersect(line[last - 1][0], line[last - 1][1], line[last][0], line[last][1])){
			last--;
			if (v[i].second[0].direction(line[last][0], line[last][1]) == 1)line[last] = v[i].second;
		}
		if (first<last)
			point[last - 1] = lineLineIntersectPoint(line[last - 1][0], line[last - 1][1], line[last][0], line[last][1]);
	}
	while (first < last && point[last - 1].direction(line[first][0], line[first][1]) == -1) last--;
	if (last - first <= 1) return 0;
	point[last] = lineLineIntersectPoint(line[first][0], line[first][1], line[last][0], line[last][1]);
	int num = unique(&*point.begin() + first, &*point.begin() + last + 1) - &point[first];
	while (num>1 && point[first] == point[first + num - 1])num--;
	memcpy(ret, &point[first], sizeof(Point)*num);
	return num;
}
int halfPlainIntersect(const Point(*p)[2], int n, Point* ret)
{
	vector<pair<double, const Point*>> v(n + 4);
	Point ext[4][2] = { { { -INF, -INF }, { INF, -INF } }, { { INF, -INF }, { INF, INF } },
	{ { INF, INF }, { -INF, INF } }, { { -INF, INF }, { -INF, -INF } } };
	for (int i = 0; i < 4; i++)
		v[i].second = ext[i];
	for (int i = 0; i < n; i++)
		v[i + 4].second = p[i];
	return halfPlainIntersectInternal(v, n + 4, ret);
}
int polygonKernel(const Point* p, int n, Point* ret)
{
	vector<pair<double, const Point*>> v;
	Point ext[2] = { p[n - 1], p[0] };
	v[0].second = ext;
	for (int i = 1; i < n; i++)
		v[i].second = &p[i - 1];
	return halfPlainIntersectInternal(v, n, ret);
}

struct NearestPointsStruct{
	Point p1, p2;
	double d2;
	vector<Point> v;
};
inline bool nearestPointsCmp(const Point& p1, const Point& p2){
	return LESS_EQUAL(p1.y, p2.y) && (LESS(p1.y, p2.y) || LESS(p1.x, p2.x));
}
void nearestPointsInternal(const Point* p, int left, int right, NearestPointsStruct& s)
{
	if (right - left < 8){
		for (int i = left; i < right; i++){
			for (int j = i + 1; i < right; j++){
				double td2 = p[j].dis2(p[i]);
				if (td2 < s.d2){
					s.d2 = td2;
					s.p1 = p[i]; s.p2 = p[j];
				}
			}
		}
		return;
	}
	int mid = (left + right) >> 1;
	nearestPointsInternal(p, left, mid, s);
	nearestPointsInternal(p, mid, right, s);
	s.v.clear();
	double l = (p[mid - 1].x + p[mid].x) / 2;
	for (int i = mid - 1; i >= left && (p[i].x - l)*(p[i].x - l) < s.d2; i++)
		s.v.push_back(p[i]);
	for (int i = mid; i<right && (p[i].x - l)*(p[i].x - l) < s.d2; i++)
		s.v.push_back(p[i]);
	sort(s.v.begin(), s.v.end(), nearestPointsCmp);
	for (unsigned int i = 0; i < s.v.size(); i++){
		for (unsigned int j = i + 1; j < s.v.size() && (p[j].y - p[i].y)*(p[j].y - p[i].y) < s.d2; j++){
			double td2 = p[j].dis2(p[i]);
			if (td2 < s.d2){
				s.d2 = td2;
				s.p1 = p[i]; s.p2 = p[j];
			}
		}
	}
}
double nearestPointsSorted(const Point* p, int n, Point& ret1, Point& ret2)
{
	NearestPointsStruct s;
	s.d2 = INF;
	s.v.reserve(n);
	nearestPointsInternal(p, 0, n, s);
	ret1 = s.p1; ret2 = s.p2;
	return sqrt(s.d2);
}
double farthestPointsConvex(const Point* p, int n, Point& ret1, Point& ret2)
{
	double d2 = 0;
	for (int i = n - 1, j = n - 2; i >0; i--){
		while (1){
			double td2 = p[i].dis2(p[j]);
			if (td2 > d2){
				d2 = td2;
				ret1 = p[i]; ret2 = p[j];
			}
			if (!j)break;
			j--;
		}
	}
	return sqrt(d2);
}
double farthestPointsSorted(const Point* p, int n, Point& ret1, Point& ret2)
{
	vector<Point> v;
	v.reserve(n);
	//convexHullSorted(p, n, &*v.begin());
	return farthestPointsConvex(&*v.begin(), v.size(), ret1, ret2);
}

int circleLineRelation(const Point& o, double r, const Point& p1, const Point& p2)
{
	double d = o.lineDis(p1, p2);
	if (LESS(d, r))return 1;
	if (LESS(r, d))return 3;
	return 2;
}
int circleCircleRelation(const Point& o1, double r1, const Point& o2, double r2)
{
	double r = o1.dis(o2);
	if (LESS(r1 + r2, r))return 4;
	if (!LESS_EQUAL(r1 + r2, r))return 3;
	double sub = abs(r1 - r2);
	if (LESS(sub, r))return 2;
	if (!LESS_EQUAL(sub, r))return 1;
	return 0;
}
bool circleLineSegIntersect(const Point& o, double r, const Point& p1, const Point& p2)
//include extream points.
{
	int t1 = p1.inCircle(o, r), t2 = p2.inCircle(o, r);
	if (t1 >= 0 || t2 >= 0)
		return t1 != 1 || t2 != 1;
	double t = o.lineRelation(p1, p2);
	if (t >= 1 || t <= 0)return false;
	Point foot = p1 + (p2 - p1)*r;
	return foot.inCircle(o, r) >= 0;
}
void circleLineIntersect(const Point& o, double r, const Point& p1, const Point& p2, Point& ret1, Point& ret2)
{
	Point foot = o.footPoint(p1, p2);
	double t = sqrt((r*r - o.dis2(foot)) / p1.dis2(p2));;
	ret1 = foot + (p2 - p1)*t;
	ret2 = foot * 2 - ret1;
}
void circleCircleIntersect(const Point& o1, double r1, const Point& o2, double r2, Point& ret1, Point& ret2)
{
	double d2 = o1.dis2(o2);
	double t1 = (r1*r1 - r2*r2) / (2 * d2) + 0.5;
	double t2 = sqrt(r1*r1 / d2 - t1*t1);
	Point foot = o1 + (o2 - o1)*t1;
	ret1 = foot + (o2 - o1).rotate90()*t2;
	ret2 = foot * 2 - ret1;
}
void Point::pointcut(const Point& o, double r, Point& ret1, Point& ret2)const
{
	double t1 = r*r / dis2(o);
	Point foot = o + (o - *this)*t1;
	double t2 = sqrt(t1 - t1*t1);
	ret1 = foot + (*this - o).rotate90()*t2;
	ret2 = foot * 2 - ret1;
}
Point Point::nearnestPoint(const Point& o, double r)const
{
	Point p = *this - o;
	double d = p.r();
	if (EQUAL(d, 0))return o;
	return o + p*(r / d);
}
//Upset the order before using this function.
double minCoveringCircle(const Point* p, int n, Point& ret)
{
	if (n == 1){ ret = p[0]; return 0; }
	double r2 = p[0].dis2(p[1]);
	ret = (p[0] + p[1]) * 0.5;
	for (int i = 2; i < n; i++){
		if (LESS(r2, ret.dis2(p[i]))){
			ret = (p[0] + p[i]) * 0.5;
			r2 = p[0].dis2(p[i]);
			for (int j = 1; j < i; j++){
				if (LESS(r2, ret.dis2(p[j]))){
					ret = (p[i] + p[j]) * 0.5;
					r2 = p[i].dis2(p[j]);
					for (int k = 0; k < j; k++){
						if (LESS(r2, ret.dis2(p[k]))){
							ret = circumcenter(p[i], p[j], p[k]);
							r2 = ret.dis2(p[k]);
						}
					}
				}
			}
		}
	}
	return sqrt(r2);
}
int unitCoveringCircle(const Point* p, int n, double r)
{
	int ret = 0;
	vector<pair<double, bool>> v;
	v.reserve(2 * n);
	double t = r*r * 4;
	for (int i = 0; i < n; i++){
		v.clear();
		int value = 0;
		for (int j = 0; j < n; j++){
			if (LESS_EQUAL(p[i].dis2(p[j]), t) && i != j){
				double a = p[j].angle(p[i]);
				double b = acos(p[i].dis(p[j]) / r / 2);
				double t1 = a - b, t2 = a + b;
				if (t1 < -PI / 2){
					if (t2 < -PI / 2){
						a += 2 * PI;
						b += 2 * PI;
					}
					else value++;
				}
				v.push_back(make_pair(t1, true));
				v.push_back(make_pair(t2, false));
			}
		}
		sort(v.begin(), v.end());
		if (value > ret)ret = value;
		for (unsigned int j = 0; j < v.size(); j++){
			if (v[j].second){
				value++;
				if (value > ret)ret = value;
			}
			else value--;
		}
	}
	return ret;
}
