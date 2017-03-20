#include "deploy.h"
#include "lib/lib_time.h"
#include <cmath>

#define MAX_NODES_NUM 1000 // 最大网络节点数量

#define MAX_GENERATION 200 // 迭代次数

int POP_SCALE = 80; // 种群规模

double Cmax; // cost -> fitness (fitness = Cmax - cost)
double Pc = 0.7; // 交叉概率
double Pm = 0.002; // 变异概率
//1.65
double c = 1.75; // 适应度变换参数

int g_is_first = 1;
/*
** return a random real in the interval
** [a, b] (also [a, b))
*/
double uniform_real(double a, double b) {
    if (g_is_first) {
        g_is_first = 0;
        srand((unsigned int)time(NULL));
    }
    return (double)rand() / ((double)RAND_MAX / (b - a)) + a;
}

/*
** return a random integer in the interval
** [a, b]
*/
int uniform_int(int a, int b) {
    if (g_is_first) {
        g_is_first = 0;
        srand((unsigned int)time(NULL));
    }

    return (int)((double)rand() / ((RAND_MAX + 1.0) / (b - a + 1.0)) + a);
}

/*
 return max element
*/
double max_ele(double a, double b) {
    if (a > b) {
        return a;
    } else {
        return b;
    }
}

/*
 return the minimum element
*/
double min_ele(double a, double b) {
    if (a < b) {
        return a;
    } else {
        return b;
    }
}

// 个体类
class Individual {

public:
    int cost; // 花费
    double fitness; // 适应度
    bitset<MAX_NODES_NUM> bitIn; // 一个个体
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

// 自定义比较函数
bool compByCost(Individual a, Individual b) {
    return a.cost < b.cost;
}

// 种群类
class Population {

public:
    vector<Individual> inVec; // 种群容器
    int bestIndex; // 当代最优个体的index
    int worstIndex; // 当代最差个体的index
    Individual bestIndividual;
    Individual worstIndividual;
    Individual everBestIndividual; // 历史最优个体
    double averageFit;
    double maxFit;
    double minFit;
    int gen; // 代数
    int iteration;
public:
    Population() { // constructor
        for (int i = 0; i < POP_SCALE; i++) {
            Individual ind = Individual(); // 初始化
            inVec.push_back(ind);
        }
    }

    // 随机选取第一代
    void generateInitalPopulation() {
        int i = 0;
        while (i < POP_SCALE) { // 每个体
            // 随机产生服务器的个数
            int serverNum = uniform_int(1, clientNum);
            //int serverNum;
            set<int> s;
            // 产生 serverNum 个随机数
            while (true) {
                int r = uniform_int(0, nodesNum-1);
                s.insert(r);
                if (s.size() == unsigned(serverNum)) {
                    break;
                }
            }
            set<int>::iterator it;
            for (it = s.begin(); it != s.end(); it++) {
                // 种群中的每个个体进行初始化
                inVec[i].bitIn.set(*it); // 置1
            }
            i++;
        }
    }

    // 计算 cost
    void calcCostValue() {
        for (int i = 0; i < POP_SCALE; i++) {
            inVec[i].cost = calCost(inVec[i].bitIn, 0);
        }
    }

    // 计算 适应度
    void calcFitnessValue() {
        double temp;
        double sumFit = 0;
        maxFit = 0;
        minFit = Cmax - double(inVec[0].cost);
        for (int i = 0; i < POP_SCALE; i++) {
            temp = Cmax - double(inVec[i].cost);
            if (temp < minFit) {
                minFit = temp;
            }
            if (temp > maxFit) {
                maxFit = temp;
            }
            sumFit += temp;
            inVec[i].fitness = temp;
        }
        // 计算当代适应度的均值
        averageFit = sumFit / (POP_SCALE * 1.0);
        /* 对适应度做线性尺度变换
            F' = a * F + b
            a = (c - 1) * F_avg / (F_max - F_avg);
            b = (F_max - c * F_avg) * F_avg / (F_max - F_avg)
        */
        double a_, b_;
        if (minFit > (c * averageFit - maxFit) / (c - 1)) {
            a_ = (c - 1) * averageFit / (maxFit - averageFit);
            b_ = ((maxFit - c * averageFit) * averageFit) / (maxFit - averageFit);
        } else {
            a_ = averageFit / (averageFit - minFit);
            b_ = minFit * averageFit / (minFit - averageFit);
        }

        for (int i = 0; i < POP_SCALE; i++) {
            temp = a_ * inVec[i].fitness + b_;
            inVec[i].fitness = temp;
        }


        /*
        // 指数
        double m = 1 + log10(MAX_GENERATION);
        cout <<"MMMMM : "<<m<<endl;
        double aa = pow(gen, m) / (averageFit + 0.001);
        cout <<"AAAAA : "<<aa<<endl;
        for (int i = 0; i < POP_SCALE; i++) {
            temp = exp(-1 * aa * inVec[i].fitness);
            cout <<"TEMP : "<<temp<<" ";
            inVec[i].fitness = temp;
        }
        cout <<endl;
        */
    }

    // find best and worst individual of this generation
    void findBestAndWorstIndividual() {

        worstIndividual.cost = inVec[0].cost;
        bestIndividual.cost = inVec[0].cost;
        worstIndex = 0;
        bestIndex = 0;
        for (unsigned int i = 1; i < inVec.size(); i++) {
            if (inVec[i].cost > worstIndividual.cost) {
                worstIndividual = inVec[i];
                worstIndex = i;
            }
            if (inVec[i].cost < bestIndividual.cost) {
                bestIndividual = inVec[i];
                bestIndex = i;
            }
        }

        // 将Cmax设为当代的最大Cost
        Cmax = double(worstIndividual.cost);

         if (gen == 1) {
            everBestIndividual = bestIndividual;
            iteration = 0;
        } else {
            if (bestIndividual.cost < everBestIndividual.cost) {
                iteration = 0;
                everBestIndividual = bestIndividual;
            } else {
                iteration++;
            }
        }
    }

    /* 精英策略
        将当代产生的最优个体与历史最优个体比较, 取最优, 替换最差
    */
    void performEvolution() {
        if (bestIndividual.cost < everBestIndividual.cost) {
            everBestIndividual = inVec[bestIndex];
        } else {
            inVec[worstIndex] = everBestIndividual;
        }
    }

    // Must Choose
    void mustChooseVec() {
        //return;
        int num = mustChoose.size();
        for (int i = 0; i < POP_SCALE; i++) {
            for (int j = 0; j < num; j++) {
                inVec[i].bitIn.set(mustChoose[j]);
            }
        }
    }

    // evaluate
    void evalutePopulation() {
        calcCostValue();
        findBestAndWorstIndividual();
        calcFitnessValue();
    }

    // 选择
    void select() {
        int i;
        double totalFitness = 0;
        // 计算总 cost, 轮盘赌
        for (i = 0; i < POP_SCALE; i++) {
            totalFitness += inVec[i].fitness;
        }

        vector<double> cumFitness;
        for (i = 0; i < POP_SCALE; i++) {
            double temp = double((inVec[i].fitness * 1.0) / (totalFitness * 1.0));
            cumFitness.push_back(temp);
        }

        // 计算累加概率
        for (i = 1; i < POP_SCALE; i++) {
            cumFitness[i] = cumFitness[i-1] + cumFitness[i];
        }

        int index;
        vector<Individual> tempVec;
        for (i = 0; i < POP_SCALE; i++) {
            //生成一个[0， 1]之间的随机浮点数
            double range = uniform_real(0, 1);
            index = 0;
            while (range > cumFitness[index]) {
                index++;
            }
            tempVec.push_back(inVec[index]);
        }
        for (int i = 0; i < POP_SCALE; i++) {
            inVec[i] = tempVec[i];
        }
    }

    // 交叉
    void crossOver() {
        int i, j;
        vector<int> index;
        int point, temp;
        double p;
        for (i = 0; i < POP_SCALE; i++) {
            index.push_back(i);
        }
        // 生成不重复随机数组
        for (i = 0; i < POP_SCALE; i++) {
            point = uniform_int(0, POP_SCALE - 1 - i);
            temp = index[i];
            index[i] = index[point + i];
            index[point + i] = temp;
        }
        bitset<BITSIZE> bita;
        bitset<BITSIZE> bitb;
        /*
        // one-point crossover
        for (i = 0; i < POP_SCALE - 1; i += 2) {
            p = uniform_real(0, 1);
            // 若小于交叉概率, 则进行交叉
            if (p < Pc) {
                point = uniform_int(0, nodesNum - 2) + 1;
                bita = inVec[index[i]].bitIn;
                bitb = inVec[index[i+1]].bitIn;
                for (j = point; j < nodesNum; j++) {
                    temp = bita[j];
                    bita[j] = bitb[j];
                    bitb[j] = temp;
                }
                // 判断服务器数量是否超过消费节点数量, 若超过, 则不交叉
                if (bita.count() <= unsigned(clientNum) && bitb.count() <= unsigned(clientNum)) {
                    inVec[index[i]].bitIn = bita;
                    inVec[index[i+1]].bitIn = bitb;
                }
            }
        }
        */

        // Uniform Crossover
        for (i = 0; i < POP_SCALE - 1; i += 2) {
            p = uniform_real(0, 1);
            if (p < Pc) {
                bita = inVec[index[i]].bitIn;
                bitb = inVec[index[i+1]].bitIn;
                for (j = 0; j < nodesNum; j++) {
                    int t = uniform_int(0, 1);
                    if (t == 1) {
                        temp = bita[j];
                        bita[j] = bitb[j];
                        bitb[j] = temp;
                    }
                }
                // 判断服务器数量是否超过消费节点数量, 若超过, 则不交叉
                if (bita.count() > 0 && bita.count() <= unsigned(clientNum) && bitb.count() > 0 && bitb.count() <= unsigned(clientNum)) {
                    inVec[index[i]].bitIn = bita;
                    inVec[index[i+1]].bitIn = bitb;
                }
            }
        }
    }

    // 变异
    void mutation() {
        double p;
        // bit mutation
        for (int i = 0; i < POP_SCALE; i++) {
            for (int j = 0; j < nodesNum; j++) {
                p = uniform_real(0, 1);
                // mutation
                if (p < Pm) {
                    // 该位取反
                    inVec[i].bitIn.flip(j);
                }
                // 变异后服务器节点数量大于消费节点数量
                if (inVec[i].bitIn.count() > unsigned(clientNum) && inVec[i].bitIn.count() > 0) {
                    inVec[i].bitIn.flip(j);
                }
            }

            /*
            // one point mutation
            int postion = uniform_int(0, nodesNum-1);
            p = uniform_real(0, 1);
            if (p < Pm) {
                inVec[i].bitIn.flip(postion);
            }
            if (inVec[i].bitIn.count() > unsigned(clientNum) && inVec[i].bitIn.count() > 0) {
                    inVec[i].bitIn.flip(postion);
            }
            */
        }
    }

    // 产生新的一代
    void generateNextPopulation() {
        select();
        crossOver();
        mutation();
    }

    // 整个进化的全过程
    void epoch() {
        //print_time("epoch begin");
        gen = 1;
        // 生成初代
        generateInitalPopulation();
        // must choose point
        mustChooseVec();
        // 计算 cost, fitness, 找出最优最差个体
        evalutePopulation();
        //show();
        //print_time("epoch end");
        while(gen < MAX_GENERATION) {
            gen++;
            // 通过选择 交叉 变异 生成下一代
            generateNextPopulation();
            mustChooseVec();
            // 评估新的一代
            evalutePopulation();
            // 精英策略
            performEvolution();
            //show();
            if (iteration > 50) {
                //break;
            }
        }
        //showBestPosition();
    }

    // 输出种群现状
    void show() {
        int i;
        int sum = 0;
        double sumFit = 0;
        double average;
        double aveFit;
        for (i = 0; i < POP_SCALE; i++) {
            sum += inVec[i].cost;
            sumFit += inVec[i].fitness;
        }
        average = double((sum * 1.0) / (POP_SCALE * 1.0));
        aveFit = double((sumFit * 1.0) / (POP_SCALE * 1.0));
        printf("Gen = %d, PopScale = %d, Average Cost = %lf, Average Fit = %lf, Min Cost = %d, Current MinCost = %d, Current worstCost = %d\n",
               gen, POP_SCALE, average, aveFit, everBestIndividual.cost, bestIndividual.cost, worstIndividual.cost);

/*
        vector<Individual>::iterator it;
        for (it = inVec.begin(); it != inVec.end(); it++) {
            (*it).show();
        }

        cout<<endl<<endl;
*/
    }

    void showBestPosition() {
        cout << "Must Choose: ";
        for (unsigned int i = 0; i < mustChoose.size(); i++) {
            cout <<mustChoose[i]<<" ";
        }
        cout <<endl;
        cout <<"Cost: "<<everBestIndividual.cost<<endl;
        cout <<"Position: ";
        for (int i = 0; i < nodesNum; i++) {
            if (everBestIndividual.bitIn[i] == 1) {
                cout <<i<<" ";
            }
        }
        cout <<endl;
    }
};
