#ifndef _INDIVIDUAL_
#define _INDIVIDUAL_
#include <bitset>
// 个体类
class Individual {

public:
    int cost; // 花费
    double fitness; // 适应度
    bitset<BITSIZE> bitIn; // 一个个体
public:
    // constructor
    Individual() {}

    // show this individual's info
    void show() {
        for (int i = 0; i < nodesNum; i++) {
            cout << bitIn[i] << " ";
        }
        cout << endl << "COST: " << cost << " FITNESS: "<< fitness << "  Server Num: " << bitIn.count() << endl;
    }
};
#endif