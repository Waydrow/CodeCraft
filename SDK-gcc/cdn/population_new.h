#include "deploy.h"
#include "lib/lib_time.h"
#include <cmath>

#define MAX_NODES_NUM 1000 // 最大网络节点数量

#define MAX_GENERATION 200 // 迭代次数

int POP_SCALE = 80; // 种群规模

double Cmax; // cost -> fitness (fitness = Cmax - cost)
double Pc = 0.7; // 交叉概率
double Pm = 0.0012; // 变异概率

double c = 1.75; // 线性变换参数
double m = 1 + log10(MAX_GENERATION); // 指数变换参数

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

bool compByFitness(Individual a, Individual b) {
    return a.fitness < b.fitness;
}

// 种群类
class Population {

public:
    vector<Individual> inVec; // 种群容器
    vector<Individual> secondVec;
    vector<Individual> buffer;
    int bestIndex; // 当代最优个体的index
    int worstIndex; // 当代最差个体的index
    Individual bestIndividual;
    Individual worstIndividual;
    Individual everBestIndividual; // 历史最优个体
    double averageFit;
    int gen; // 代数
    int iteration;
public:
    Population() { // constructor
        for (int i = 0; i < POP_SCALE; i++) {
            Individual ind = Individual(); // 初始化
            inVec.push_back(ind);
        }
    }

    void generateIndividual(Individual &in) {
        // 随机产生服务器的个数
        int serverNum = uniform_int(1, clientNum - int(mustChoose.size()));
        //int serverNum;
        set<int> s;
        for (unsigned int j = 0; j < mustChoose.size(); j++) {
            s.insert(mustChoose[j]);
        }
        // 产生 serverNum 个随机数
        while (true) {
            int r = uniform_int(0, nodesNum-1);
            s.insert(r);
            if (int(s.size()) == serverNum + int(mustChoose.size())) {
                break;
            }
        }
        set<int>::iterator it;
        in.bitIn.reset();
        for (it = s.begin(); it != s.end(); it++) {
            in.bitIn.set(*it);
        }
    }

    // 随机选取第一代
    void generateInitalPopulation() {
        int i = 0;
        while (i < POP_SCALE) { // 每个体
            generateIndividual(inVec[i]);
            if (!meetGH(inVec[i], i, inVec)) {
                continue;
            }
            i++;
        }
    }

    // 检查该个体与整个种群的海明距离
    bool meetGH(Individual a, unsigned int index, vector<Individual> vec) {
        for (unsigned int i = 0; i < vec.size(); i++) {
            if (index == i) {
                continue;
            }
            int GH = (a.bitIn ^ vec[i].bitIn).count();
            if (GH < 3) {
                return false;
            }
        }
        return true;
    }

    // 计算 cost
    void calcCostValue(vector<Individual> &vec) {
        double sumFit = 0;
        for (unsigned int i = 0; i < vec.size(); i++) {
            vec[i].cost = calCost(vec[i].bitIn, 0, true);
            sumFit += double(vec[i].cost);
        }
        averageFit = sumFit / (POP_SCALE * 1.0);
    }

    // 计算 适应度
    void calcFitnessValue(vector<Individual> &vec) {

        // 指数
        double aa = pow(gen, 1 / m) / (averageFit + 0.001);
        for (unsigned int i = 0; i < vec.size(); i++) {
            vec[i].fitness = exp(-1 * aa * double(vec[i].cost));
        }
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


    // Must Choose
    void mustChooseVec(vector<Individual> &vec) {
        //return;
        int num = mustChoose.size();
        for (unsigned int i = 0; i < vec.size(); i++) {
            for (int j = 0; j < num; j++) {
                vec[i].bitIn.set(mustChoose[j]);
            }
        }
    }

    void sortByFitness(vector<Individual> &vec) {
        // 由小到大排列
        sort(vec.begin(), vec.end(), compByFitness);
    }

    // evaluate
    void evalutePopulation(vector<Individual> &vec) {
        calcCostValue(vec);
        //findBestAndWorstIndividual();
        calcFitnessValue(vec);
        sortByFitness(vec);

        if (gen == 1) {
            everBestIndividual = vec[vec.size() - 1];
        } else {
            if (vec[vec.size() - 1].cost < everBestIndividual.cost) {
                everBestIndividual = vec[vec.size() - 1];
            }
        }
    }

    // 选择
    vector<Individual> select(vector<Individual> vec) {

        // 轮盘赌
        int i;
        double totalFitness = 0;
        // 计算总 cost, 轮盘赌
        for (i = 0; i < POP_SCALE; i++) {
            totalFitness += vec[i].fitness;
        }

        vector<double> cumFitness;
        for (i = 0; i < POP_SCALE; i++) {
            double temp = double((vec[i].fitness * 1.0) / (totalFitness * 1.0));
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
            tempVec.push_back(vec[index]);
        }

        return tempVec;

    }

    // 交叉
    void crossOverAndMutation(vector<Individual> &vec, double averageGH, double maxFit, double avgFit) {
        unsigned int i, j;
        vector<int> index;
        int point, temp;
        double p;
        for (i = 0; i < vec.size(); i++) {
            index.push_back(i);
        }
        // 生成不重复随机数组
        for (i = 0; i < vec.size(); i++) {
            point = uniform_int(0, vec.size() - 1 - i);
            temp = index[i];
            index[i] = index[point + i];
            index[point + i] = temp;
        }
        bitset<BITSIZE> bita;
        bitset<BITSIZE> bitb;


        // Uniform Crossover
        for (i = 0; i < vec.size() - 1; i += 2) {
            double alpha = 0.3;
            int index1 = index[i];
            int index2 = index[i+1];
            int GH = (vec[index1].bitIn ^ vec[index2].bitIn).count();
            double t_ = alpha + (GH - averageGH) / (0.001 + averageGH);
            if (t_ < 0) {
                Pc = 0;
            } else if (t_ <= 1) {
                Pc = t_;
            } else {
                Pc = 1;
            }
            p = uniform_real(0, 1);
            if (p < Pc) {
                bita = vec[index1].bitIn;
                bitb = vec[index2].bitIn;
                for (j = 0; j < unsigned(nodesNum); j++) {
                    int t = uniform_int(0, 1);
                    if (t == 1) {
                        temp = bita[j];
                        bita[j] = bitb[j];
                        bitb[j] = temp;
                    }
                }

                // 判断服务器数量是否超过消费节点数量, 若超过, 则不交叉
                if (bita.count() > 0 && bita.count() <= unsigned(clientNum) && bitb.count() > 0 && bitb.count() <= unsigned(clientNum)) {
                    vec[index1].bitIn = bita;
                    vec[index2].bitIn = bitb;
                }
            }
            mutation(vec[index1], GH, averageGH, maxFit, avgFit);
            mutation(vec[index2], GH, averageGH, maxFit, avgFit);
        }
    }

    // 变异
    void mutation(Individual &in, int GH, int averageGH, double maxFit, double avgFit) {
        double p;
        double beta = 0.005;
        Pm = beta / ((maxFit - avgFit) * averageGH + (GH - averageGH) + 0.001);
        //cout << "PM: "<<Pm<<endl;
        // bit mutation
        for (int j = 0; j < nodesNum; j++) {
            p = uniform_real(0, 1);
            // mutation
            if (p < Pm) {
                // 该位取反
                in.bitIn.flip(j);
            }
            // 变异后服务器节点数量大于消费节点数量
            if (in.bitIn.count() > unsigned(clientNum) || in.bitIn.count() == 0) {
                in.bitIn.flip(j);
            }
        }
    }

    // 产生种群2, 取前 0.1*N 个进入种群 2
    void copyToSecond() {
        secondVec.clear();
        int num = int(0.1 * POP_SCALE);
        for (int i = POP_SCALE - 1; i >= POP_SCALE - num; i--) {
            secondVec.push_back(inVec[i]);
            inVec.pop_back();
        }
    }

    // 将种群 1 中剩余的 0.9*N 个体复制到缓冲区
    void copyToBuffer() {
        buffer.clear();
        double sumFit = 0;
        for (unsigned int i = 0; i < inVec.size(); i++) {
            sumFit += inVec[i].fitness;
        }
        int num = int(0.9 * POP_SCALE);
        for (int i = 0; i < num; i++) {
            int nextCount = int (inVec[i].fitness / sumFit * num * 1.0 + 0.5);
            for (int j = 0; j < nextCount; j++) {
                buffer.push_back(inVec[i]);
            }
        }

        /*
        while (buffer.size() < unsigned(num)) {
            Individual in = Individual();
            generateIndividual(in);
            buffer.push_back(in);
        }
        */

    }

    void generateSecondPopulation() {
        copyToSecond();
        copyToBuffer();
        int maxGH;
        double averageGH;
        while (true) {
            calGH(buffer, averageGH, maxGH);
            if (averageGH < 0.01) {
                Individual in = Individual();
                generateIndividual(in);
                int position = uniform_int(0, buffer.size() - 1);
                buffer[position] = in;
            } else {
                break;
            }
        }
        double sumFit = 0;
        double maxFit = 0;
        double avgFit;
        for (unsigned int i = 0; i < buffer.size(); i++) {
            double temp = buffer[i].fitness;
            if (temp > maxFit) {
                maxFit = temp;
            }
            sumFit += temp;
        }
        avgFit = sumFit / double(buffer.size());
        crossOverAndMutation(buffer, averageGH, maxFit, avgFit);
        for (unsigned int i = 0; i < buffer.size(); i++) {
            secondVec.push_back(buffer[i]);
        }

    }

    void calGH(vector<Individual> &vec, double &averageGH, int &maxGH) {
        int vecSize  = int(vec.size());
        int sumGH = 0;
        maxGH = 0;
        int num = 0;
        for (int i = 0; i < vecSize - 1; i++) {
            for (int j = i + 1; j < vecSize; j++) {
                int temp = (vec[i].bitIn ^ vec[j].bitIn).count();
                if (temp > maxGH) {
                    maxGH = temp;
                }
                sumGH += temp;
                num++;
            }
        }
        averageGH = (sumGH * 1.0) / (num * 1.0);
    }

    // 整个进化的全过程
    void epoch() {
        gen = 1;
        generateInitalPopulation(); // 产生第一代
        evalutePopulation(inVec); // 评价第一代并排序
        show(inVec);
        while(gen < MAX_GENERATION) {

            generateSecondPopulation(); // 产生第二代
            mustChooseVec(secondVec);
            evalutePopulation(secondVec);
            inVec = secondVec;
            gen++;
            show(inVec);
        }
        showBest();
    }

    void showBest() {
        cout<<endl<<endl<<"BEST: "<<everBestIndividual.cost<<endl;
    }

    // 输出种群现状
    void show(vector<Individual> vec) {
        int i;
        int sum = 0;
        double sumFit = 0;
        double average;
        double aveFit;
        for (i = 0; i < POP_SCALE; i++) {
            sum += vec[i].cost;
            sumFit += vec[i].fitness;
        }
        average = double((sum * 1.0) / (POP_SCALE * 1.0));
        aveFit = double((sumFit * 1.0) / (POP_SCALE * 1.0));
        printf("Gen = %d, PopScale = %d, Average Cost = %lf, Average Fit = %lf, Min Cost = %d, Current MinCost = %d, Current worstCost = %d\n",
               gen, POP_SCALE, average, aveFit, everBestIndividual.cost, vec[vec.size() - 1].cost, vec[0].cost);

/*
        vector<Individual>::iterator it;
        for (it = vec.begin(); it != vec.end(); it++) {
            (*it).show();
        }

        cout<<endl<<endl;
*/
    }
};
