#include "deploy.h"
#include <stdio.h>
#include <iostream>
#include <string>



using namespace std;

int nodesNum; // 网络节点数量
int linkNum; // 网络链路数量
int clientNum; // 消费节点数量
int deployCost; // 服务器部署成本


string relStr;

//你要完成的功能总入口
/* topo[MAX_EDGE_NUM] 中每一项存储输入文件中的一行
   line_num 为输入文件中的总行数
   filename 为输出文件名
*/

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename) {

    sscanf(topo[0], "%d%d%d", &nodesNum, &linkNum, &clientNum);
    sscanf(topo[2], "%d", &deployCost);

    char tempSSS[100];
    sprintf(tempSSS, "%d\n\n", clientNum);
    relStr += tempSSS;
    int temp = linkNum + 5;
    for (int i = temp; i < temp + clientNum; i++) {
        int a, b, c;
        sscanf(topo[i], "%d%d%d", &a, &b, &c);
        char ttts[200];
        sprintf(ttts, "%d %d %d", b, a, c);
        relStr += ttts;
        if (i < temp + clientNum - 1) {
            relStr += "\n";
        }
    }

    // 直接调用输出文件的方法输出到指定文件中(ps请注意格式的正确性，如果有解，第一行只有一个数据；第二行为空；第三行开始才是具体的数据，数据之间用一个空格分隔开)
    write_result(relStr.c_str(), filename);
}
