template<class T>
struct LeftistTree {
    struct Data {
        T key;
        int l, r, dist;
    } D[MX << 1];
    int rear, root;

    void init() {
        rear = root = 0;
        D[0].dist = -1;
    }
    LeftistTree() {
        init();
    }

    int New(T key) {
        rear++;
        D[rear].l = D[rear].r = 0;
        D[rear].key = key;
        D[rear].dist = 0;
        return rear;
    }
    int merge(int r1, int r2) {
        if(!r1) return r2;
        if(!r2) return r1;
        if(D[r2].key < D[r1].key) {
            swap(r1, r2);
        }
        D[r1].r = merge(D[r1].r, r2);
        if(D[D[r1].l].dist < D[D[r1].r].dist) {
            swap(D[r1].l, D[r1].r);
        }
        D[r1].dist = D[D[r1].r].dist + 1;
        return r1;
    }
    void push(int &rt, T key) {
        rt = merge(rt, New(key));
    }
    T pop(int &rt) {
        T ret = D[rt].key;
        rt = merge(D[rt].l, D[rt].r);
        return ret;
    };
};
