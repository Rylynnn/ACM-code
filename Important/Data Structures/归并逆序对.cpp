while(l<=mid && r<=y)
    {
        if(a[l]>a[r])
        {
            temp[i++]=a[r++];
            ans+=mid-l+1;
            continue;
        }
        temp[i++]=a[l++];
    }
