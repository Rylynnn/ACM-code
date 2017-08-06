int Gao(char a[],int len) {
  int i = 0,j = 1,k = 0;
  while (i < len && j < len && k < len) {
    int cmp = a[(j+k)%len]-a[(i+k)%len];
    if (cmp == 0)
      k++;
    else {
      if (cmp > 0)
        j += k+1;
      else
        i += k+1;
      if (i == j)	j++;
      k = 0;
    }
  }
  return min(i,j);
}

int MinimumRepresentation(int *s, int l)    
{    
    int i,j,k;  
    i=0;j=1;k=0;  
    while(i<l&&j<l)  
    {  
        k=0;  
        while(s[i+k]==s[j+k]&&k<l) k++;  
        if(k==l) return i;  
        if(s[i+k]>s[j+k])   
         if(i+k+1>j) i=i+k+1;  
         else i=j+1;  
        else if(j+k+1>i) j=j+k+1;  
        else  j=i+1;   
    }  
    if(i<l) return i;  
    else return j;  
}   
