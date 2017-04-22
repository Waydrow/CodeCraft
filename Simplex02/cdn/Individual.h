#ifndef _Individual
#define _Individual

// 个体类
class Individual {

public:
    int cost; // 花费
    double fitness; // 适应度
    //bitset<MAX_NODES_NUM> bitIn; // 一个个体
    int gen[BITSIZE];
public:
    // constructor
    Individual() {
        memset(gen, 0, sizeof(gen));
    }

    ~Individual() {}

    int count() {
        int sum = 0;
        for (int i=0; i<nodesNum; i++) {
            if (gen[i]) sum++;
        }
        return sum;
    }

    friend Individual operator |(const Individual &a,const Individual &b){
        Individual c;
        for(int i=0;i<nodesNum;i++){
            c.gen[i]=max(a.gen[i],b.gen[i]);
        }
        return c;
    }

    // show this individual's info
    void show() {
        for (int i = 0; i < nodesNum; i++) {
            cout << gen[i] << " ";
        }
        cout << endl << "COST: " << cost << " FITNESS: "<< fitness << "  Server Num: " << count() << endl;
    }
};

#endif