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
#include <queue>


//..xt:增广路跑费用流[SPFA]（可heap优化dij优化）
//..xt:每次需要重新建图 （可memcpy优化）
//..xt:需要保证源点编号最小，汇点编号最大

//..xt:定义最大点数边数
#define MAXN 2000
#define MAXM (52000*8)
#define INF 10000000
#define BITSIZE 1000


using namespace std;

int nodesNum; // 网络节点数量
int linkNum; // 网络链路数量
int clientNum; // 消费节点数量
int deployCost; // 服务器部署成本

clock_t start;

string relStr;

#include "normal.h"
#include "dinic_cost.h"
#include "population.h"



//你要完成的功能总入口
/* topo[MAX_EDGE_NUM] 中每一项存储输入文件中的一行
   line_num 为输入文件中的总行数
   filename 为输出文件名
*/

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename) {
    start = clock();

    sscanf(topo[0], "%d%d%d", &nodesNum, &linkNum, &clientNum);
    sscanf(topo[2], "%d", &deployCost);
    readData(topo,nodesNum,linkNum,clientNum);//将数据从缓存中读到数组中
    getMustChoose();//挑选出必须选择的直接相连的点

/*
    int maxTimes = 0;
    double tpc, tpm;
    int cl = 0;
    for (double i = 0.5; i < 0.9; i+=0.1) {
        Pc = i;
        for (double j = 0.001; j < 0.009; j+=0.001) {
            Pm = j;
            int temp = 0;
            for (int k = 0; k < 20; k++) {
                Population p = Population();
                p.epoch();
                mp.clear();
                if (p.everBestIndividual.cost == 2111) {
                    temp++;
                }
            }
            printf("%-4d  %-4d   %.2lf, %.4lf\n", cl, temp, i, j);
            cl++;
            if (temp > maxTimes) {
                maxTimes = temp;
                tpc = i;
                tpm = j;
            }
        }
    }
    cout <<endl;
    cout <<"BEST"<<endl;
    printf("%d  %d   %.2lf, %.4lf\n", cl, maxTimes, tpc, tpm);
*/
/*
        for (int i = 0; i < 20; i++) {
            Population p = Population();
            p.epoch();
            mp.clear();
            cout <<"Min Cost: "<<p.everBestIndividual.cost<<endl;
        }
*/

    Population p = Population();//初始化种群参数
    p.epoch();//整个迭代环境
    //printf("%.5lf\n",(double)mytime/mynum/CLOCKS_PER_SEC);
    bitset<BITSIZE>rel = p.everBestIndividual.bitIn;//取得最优个体DNA
    calCost(rel,1,false);//目的并非计算cost，而是构造网络环境，从而计算有哪些clients未满足
    rel|=getBetter();//将之前relDNA中未布置服务器的点部署服务器，得到真实DNA
    printf("%d\n",calCost(rel,1,false));//对真实DNA，构造网络环境

    cout <<"Server Num: "<<rel.count()<<endl;
    printRel(rel);//打印结果

    // 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
    write_result(relStr.c_str(), filename);
    //printf("%.5lf\n",zong/jishu);
}
