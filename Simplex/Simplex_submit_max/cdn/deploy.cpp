#include "deploy.h"
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <queue>
#include <bitset>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <set>
#include <algorithm>
#include <map>
#include <utility>
#include <cmath>
#include <fstream>
enum dataProperty :int{
    num_line=0,
    servers_line,
    nodes_line,
    links_line,
    consumtionNodes_line
};

//..xt:增广路跑费用流[SPFA]（可heap优化dij优化）
//..xt:每次需要重新建图 （可memcpy优化）
//..xt:需要保证源点编号最小，汇点编号最大

//..xt:定义最大点数边数
#define MAXN 2000
#define MAXM (52000*8)
#define INF 10000000
#define BITSIZE 1300


using namespace std;

int nodesNum; // 网络节点数量
int linkNum; // 网络链路数量
int clientNum; // 消费节点数量
int deployCost; // 服务器部署成本
int serverLevel; // 服务器档次数量

clock_t start;

string relStr;
struct edge {
    edge *next,*op;
    int t,c,v,bf;
} ES[MAXM],*V[MAXN],*cur[MAXN];
int myEdgeNum[MAXN];


// 全局时间种子标志
int g_srand = 1;

const bool * firstMap;

// get a real random number in [a, b]
inline double get_random_real(double a, double b) {
    if (g_srand) {
        g_srand = 0;
        srand((unsigned int)time(NULL));
    }
    return (double)rand() / ((double)RAND_MAX / (b - a)) + a;
}

// get a int random number in [a, b]
inline int get_random_int(int a, int b) {
    if (g_srand) {
        g_srand = 0;
        srand((unsigned int)time(NULL));
    }

    return (int)((double)rand() / ((RAND_MAX + 1.0) / (b - a + 1.0)) + a);
}

#include "normal.h"
#include "dinic_cost.h"
#include "MinCostFlowSolution.h"
#include "individual.h"
#include "SA.h"
#include "population.h"



//你要完成的功能总入口
/* topo[MAX_EDGE_NUM] 中每一项存储输入文件中的一行
   line_num 为输入文件中的总行数
   filename 为输出文件名
*/

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename) {
    start = clock();

    sscanf(topo[0], "%d%d%d", &nodesNum, &linkNum, &clientNum);
    MinCostFlowSolution *mcf = new MinCostFlowSolution();
    mcf->ReadData(topo, nodesNum, linkNum, clientNum);
    firstMap = mcf->getFirstMap();
    readData(topo, nodesNum, linkNum, clientNum);//将数据从缓存中读到数组中
    initNodesP();
    buildBasicGraph();
    //getMustChoose();//挑选出必须选择的直接相连的点
    Population p = Population(mcf);
    p.epoch();
    bitset<BITSIZE>rel = p.everBestIndividual.bitIn; //取得最优个体DNA
    // bitset<BITSIZE> rel = SA(p.everBestIndividual.bitIn, mcf);
    //calCost(rel,1,false); //目的并非计算cost，而是构造网络环境，从而计算有哪些clients未满足
    //rel|=getBetter(); //将之前relDNA中未布置服务器的点部署服务器，得到真实DNA
    vector<int>vrel = mcf->getLastSolution(rel);
    printf("Cost: %d\n",calCost(vrel,1,false));//对真实DNA，构造网络环境

    //cout <<"Server Num: "<<rel.count()<<endl;

    printRel(vrel);//打印结果
    delete mcf;
    // 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
    write_result(relStr.c_str(), filename);
}
