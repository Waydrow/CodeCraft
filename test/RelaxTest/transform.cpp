#include<bits/stdc++.h>
#include<time.h>
#include<fstream>
#define MAXN 2000
#define MAXM (52000*8)
#define INF 10000000
#define BITSIZE 1000
using namespace std;
char topo[10000][100];
bitset<BITSIZE>gene;
ofstream out;
int nodesNum,linkNum,clientNum,deployCost;
int sumCost1;
#include "normal.h"
#include "dinic_cost.h"
void getTopo(){
    int i=0;
    while(fgets(topo[i++],100,stdin));
    printf("sucess %d lines\n",i);
}
void readData() {
    sscanf(topo[0], "%d%d%d", &nodesNum, &linkNum, &clientNum);
    sscanf(topo[2], "%d", &deployCost);
    int lineAdd=4,numAdd=1;
    int st,ed,bandW,cost;
    for(int i=0; i<linkNum; i++) {
        sscanf(topo[i+lineAdd],"%d%d%d%d",&st,&ed,&bandW,&cost);
        links[i][0]=st;
        links[i][1]=ed;
        links[i][2]=bandW;
        links[i][3]=cost;
    }
    lineAdd=linkNum+5;
    numAdd+=nodesNum;
    int n_client,n_node,need;
    for(int i=0; i<clientNum; i++) {
        sscanf(topo[i+lineAdd],"%d%d%d",&n_client,&n_node,&need);
        myNodes[n_client]=n_node;
        consumptionNodes[i][0]=n_client;
        consumptionNodes[i][1]=n_node;
        consumptionNodes[i][2]=need;
        sumCost1+=need;
    }
}
int main()
{
    out.open("_case0.txt");
    getTopo();
    readData();

    gene.reset();
    gene.set(2);
    gene.set(5);

    for(int i=20;i<=59;i++)gene.set(i);
    for(int i=127;i<=156;i++)gene.set(i);
    for(int i=200;i<=230;i++)gene.set(i);

    int jiedian=nodesNum+clientNum+1;
    int bian=linkNum*2+gene.count()+clientNum*2;
    out << "p min "<<jiedian<<" "<<bian<<endl;
    int sst=clock();
    int www=calCost(gene,true,false);
    int eed=clock();
    printf("time:%.5lf\n",(double)(eed-sst)/CLOCKS_PER_SEC);
    printf("myCost=%d\n",www);
    out << "n "<<1<<" "<<sumCost1<<endl;
    for(int i=2;i<2+nodesNum;i++){
        out<<"n "<<i<<" 0"<<endl;
    }
    for(int i=2+nodesNum;i<2+nodesNum+clientNum;i++){
        out<<"n "<<i<<" "<<-consumptionNodes[i-2-nodesNum][2]<<endl;
    }
    for(int i=nodesNum+1;i<=nodesNum+clientNum;i++){
        out<<"a 1 "<<i+1<<" "<<0<<" "<<consumptionNodes[i-nodesNum-1][2]<<" "<<400000<<endl;
        //out<<"a 1 "<<i+1<<" "<<0<<" "<<4000<<" "<<400000<<endl;
    }
}
