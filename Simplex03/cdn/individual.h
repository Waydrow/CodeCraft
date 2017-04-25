#ifndef _INDIVIDUAL_
#define _INDIVIDUAL_
#include <bitset>
// 个体类
class Individual {

public:
    int cost; // 花费
    double fitness; // 适应度
    bitset<BITSIZE> bitIn; // 一个个体
    int gen[1300];
public:
    // constructor
    Individual() {}

    // 由 0/1 串置换为十进制编码
    void updateRealPlan() {
        int nn = 0;
        for (int i = 0; i < nodesNum * 3; i += 3) {
            int x = 0;
            for (int j = 2; j >= 0; j--) {
                x += bitIn[i+j] * pow(2, 2-j);
            }
            gen[nn++] = x;
        }
    }

    // 由十进制编码置换为 0/1 串
    void updateBinaryPlan() {
        int nn = 0;
        for (int i = 0; i < nodesNum; i++) {
            int x = gen[i];
            for (int j = 2; j >= 0; j--) {
                if (x & (1 << j)) {
                    bitIn[nn++] = 1;
                } else {
                    bitIn[nn++] = 0;
                }
            }
        }
    }

    // return server num
    int count() {
        int x = 0;
        for (int i = 0; i < nodesNum; i++) {
            if (gen[i]) {
                x++;
            }
        }
        return x;
    }

    // show this individual's info
    void show() {
        for (int i = 0; i < nodesNum; i++) {
            cout << bitIn[i] << " ";
        }
        cout << endl << "COST: " << cost << " FITNESS: "<< fitness << "  Server Num: " << bitIn.count() << endl;
    }
};
#endif