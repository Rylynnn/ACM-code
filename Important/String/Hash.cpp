#define max_strlen 15
#define maxn 149993
struct hash_point
{
	char origin[max_strlen],dat[max_strlen];
	hash_point* next;
}hash_table[maxn],*hash_head[maxn];
int tot;
int get_hash_value(char *key)  
{  
    unsigned int h=0;  
    while(*key)  
    {  
        h=(h<<4)+(*key++);
        unsigned int g=h&0xf0000000L;  
        if(g) h^=g>>24;  
        h&=~g;  
    }  
    return h%maxn;
}
void add_hash_point(char* s,char *v)
{
	tot++;
	strcpy(hash_table[tot].dat,v);
	strcpy(hash_table[tot].origin,s);
    int ths=get_hash_value(s);
    hash_table[tot].next=hash_head[ths];
    hash_head[ths]=hash_table+tot;
}
