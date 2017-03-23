#include<bits/stdc++.h>
using namespace std;
int main()
{
    int nodesNum,linkNum,clientNum,cost,liu=10000,dc=5;
    nodesNum=63;
    linkNum=62;
    clientNum=32;
    cost=505;
    printf("%d %d %d\n\n",nodesNum,linkNum,clientNum);
    printf("%d\n\n",cost);
    int num=0;
    for(int i=0;i<64;i++){
        if(i*2+1>62)continue;
        //num++;
        printf("%d %d %d %d\n",i,i*2+1,liu,dc);
        if((i+1)*2>62)continue;
        //num++;
        printf("%d %d %d %d\n",i,i*2+2,liu,dc);
    }
    puts("");
    //printf("%d\n",num);
    for(int i=0;i<clientNum;i++){
        printf("%d %d 10\n",i,i+31);
    }
}
