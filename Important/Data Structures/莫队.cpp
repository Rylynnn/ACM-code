struct Query{
	int pos[2], i;
	bool operator < (const Query& q)const{
		int t1 = pos[0] / block, t2 = q.pos[0] / block;
		return t1 < t2 || (t1 == t2 && pos[1] < q.pos[1]);
	}
}q[200001];

