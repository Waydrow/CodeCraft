#include "deploy.h"
#include "lib/lib_time.h"
#include <cmath>

#define MAX_NODES_NUM 1000 // 最大网络节点数量

int MAX_GENERATION; // 迭代次数

double ALL_TIME = 86; // 程序最长运行时间限制 单位 s

int POP_SCALE; // 种群规模

int N_POP = 100; // 记忆前 N 个当前种群个体

double L = 3.2; // 最小 Hamming 距离

double Pc = 0.7; // 交叉概率
// 0.0012
double Pm = 0.0012; // 变异概率

double c = 1.75; // 线性变换参数
double m; // 指数变换参数

// 全局时间种子标志
int g_srand = 1;

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

// 个体类
class Individual {

public:
    int cost; // 花费
    double fitness; // 适应度
    bitset<MAX_NODES_NUM> bitIn; // 一个个体
public:
    // constructor
    Individual() {

    }

    // show this individual's info
    void show() {
        for (int i = 0; i < nodesNum; i++) {
            cout << bitIn[i] << " ";
        }
        cout << endl << "COST: " << cost << " FITNESS: "<< fitness << "  Server Num: " << bitIn.count() << endl;
    }
};

// 自定义比较函数
inline bool compByCost(Individual a, Individual b) {
    return a.cost < b.cost;
}

inline bool compByFitness(Individual a, Individual b) {
    return a.fitness > b.fitness;
}

// 种群类
class Population {

public:
    vector<Individual> inVec; // 种群容器
    vector<Individual> parentVec;
    vector<Individual> nVec;
    int bestIndex; // 当代最优个体的index
    int worstIndex; // 当代最差个体的index
    Individual bestIndividual;
    Individual worstIndividual;
    Individual everBestIndividual; // 历史最优个体
    double averageFit;
    double minFit;
    double maxFit;
    int gen; // 代数
    int iteration;
    MinCostFlowSolution *mcf;

    ~Population() {
        // cout << "123" << endl;
    }
public:
    Population(MinCostFlowSolution *mcfs) { // constructor
        mcf = mcfs;

        if (nodesNum < 200) { // 对于初级 case
            POP_SCALE = 200;
            MAX_GENERATION = 5;
            Pm = 0.0012;
        } else if (nodesNum < 500) { // 对于中级 case
            POP_SCALE = 200;
            MAX_GENERATION = 5;
            Pm = 0.0005;
            N_POP = 100;
            L = 1.84;
        } else { // 对于高级 case
            POP_SCALE = 200;
            MAX_GENERATION = 300;
            Pm = 0.0002;
            L = 2.5;
        }
        // 设置指数变换参数
        m = 1 + log10(MAX_GENERATION);
        //printf("%d %d %lf\n", POP_SCALE, MAX_GENERATION, m);


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
            int serverNum = get_random_int(1, clientNum - int(mustChoose.size()));
            //int serverNum;
            set<int> s;
            for (unsigned int j = 0; j < mustChoose.size(); j++) {
                s.insert(mustChoose[j]);
            }
            // 产生 serverNum 个随机数
            while (true) {
                int r = get_random_int(0, nodesNum-1);
                //if (notChoose.find(r) != notChoose.end()) continue;
                s.insert(r);
                if (int(s.size()) == serverNum + int(mustChoose.size())) {
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

    // 计算 cost, 与适应度指数变换对应
    void calcCostValueForExp() {
        double sumFit = 0;
        for (int i = 0; i < POP_SCALE; i++) {
            bool isThreeLevel = false;
            if (nodesNum > 500 && gen < 60) {
                isThreeLevel = true;
            }
            pair<pair<long long,long long>, bitset<BITSIZE> > pair1 = mcf->CalCost(inVec[i].bitIn, 0, isThreeLevel);
            inVec[i].cost = pair1.first.first;

            sumFit += double(inVec[i].cost);

            if (inVec[i].cost < everBestIndividual.cost) {
                everBestIndividual = inVec[i];
            }
            inVec[i].bitIn = pair1.second;
            inVec[i].cost = pair1.first.second;
            /*
            if (gen <= 60 && nodesNum > 500){
                inVec[i].bitIn = pair1.second;
                inVec[i].cost = pair1.first.second;
            }*/
            /*
            if (gen <= 30 && nodesNum > 200 && nodesNum < 500){
                inVec[i].bitIn = pair1.second;
                inVec[i].cost = pair1.first.second;
            }*/
        }
        averageFit = sumFit / (POP_SCALE * 1.0);
    }


    // 计算适应度, 指数变换
    void calcFitnessValueExp() {
        // 指数
        double aa = pow(gen, 1 / m) / (averageFit + 0.001);
        for (int i = 0; i < POP_SCALE; i++) {
            inVec[i].fitness = exp(-1 * aa * double(inVec[i].cost));
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

    // 精英策略: 替换最差个体
    inline void performEvolution() {
        inVec[worstIndex] = everBestIndividual;
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
        /*
        for (int i = 0; i < POP_SCALE; i++) {
            set<int>::iterator it;
            for (it = notChoose.begin(); it != notChoose.end(); it++) {
                inVec[i].bitIn.reset(*it);
            }
        }
        */
    }


    // evaluate
    void evalutePopulation() {
        calcCostValueForExp();
        //calcCostValueForLinear();
        findBestAndWorstIndividual();
        calcFitnessValueExp();
    }

    // 轮盘赌
    void selectByRoulette() {

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
            double range = get_random_real(0, 1);
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

    // 随机联赛选择
    void selectByLeague() {

        parentVec = inVec;
        inVec.clear();
        int N = 2;
        while (inVec.size() != unsigned(POP_SCALE)) {
            for (int i = 0; i < N; i++) {
                int position1 = get_random_int(0, POP_SCALE-1);
                int position2 = get_random_int(0, POP_SCALE-1);
                if (parentVec[position1].fitness > parentVec[position2].fitness) {
                    inVec.push_back(parentVec[position1]);
                } else {
                    inVec.push_back(parentVec[position2]);
                }
            }
        }
    }

    // 均匀交叉
    void crossOverUniform() {
        int i, j;
        vector<int> index;
        int point, temp;
        double p;
        for (i = 0; i < POP_SCALE; i++) {
            index.push_back(i);
        }
        // 生成不重复随机数组
        for (i = 0; i < POP_SCALE; i++) {
            point = get_random_int(0, POP_SCALE - 1 - i);
            temp = index[i];
            index[i] = index[point + i];
            index[point + i] = temp;
        }
        bitset<BITSIZE> bita;
        bitset<BITSIZE> bitb;

        // Uniform Crossover
        for (i = 0; i < POP_SCALE - 1; i += 2) {
            p = get_random_real(0, 1);
            if (p < Pc) {
                bita = inVec[index[i]].bitIn;
                bitb = inVec[index[i+1]].bitIn;
                for (j = 0; j < nodesNum; j++) {
                    int t = get_random_int(0, 1);
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

    // 单点交叉
    void crossOverOnePoint() {
        int i, j;
        vector<int> index;
        int point, temp;
        double p;
        for (i = 0; i < POP_SCALE; i++) {
            index.push_back(i);
        }
        // 生成不重复随机数组
        for (i = 0; i < POP_SCALE; i++) {
            point = get_random_int(0, POP_SCALE - 1 - i);
            temp = index[i];
            index[i] = index[point + i];
            index[point + i] = temp;
        }
        bitset<BITSIZE> bita;
        bitset<BITSIZE> bitb;

        for (i = 0; i < POP_SCALE - 1; i += 2) {
            p = get_random_real(0, 1);
            // 若小于交叉概率, 则进行交叉
            if (p < Pc) {
                point = get_random_int(0, nodesNum - 2) + 1;
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
    }

    // 均匀变异
    void mutationUniform() {
        double p;
        // bit mutation
        for (int i = 0; i < POP_SCALE; i++) {
            for (int j = 0; j < nodesNum; j++) {
                p = get_random_real(0, 1);
                // mutation
                if (p < Pm) {
                    // 该位取反
                    inVec[i].bitIn.flip(j);
                }
                // 变异后服务器节点数量大于消费节点数量
                if (inVec[i].bitIn.count() > unsigned(clientNum) || inVec[i].bitIn.count() == 0) {
                    inVec[i].bitIn.flip(j);
                }
            }
        }
    }

    // 单点变异
    void mutationOnePoint() {
        double p;
        // bit mutation
        for (int i = 0; i < POP_SCALE; i++) {
            int postion = get_random_int(0, nodesNum-1);
            p = get_random_real(0, 1);
            if (p < Pm) {
                inVec[i].bitIn.flip(postion);
            }
            if (inVec[i].bitIn.count() > unsigned(clientNum) && inVec[i].bitIn.count() > 0) {
                inVec[i].bitIn.flip(postion);
            }
        }
    }

    // 记忆当前种群的前 N_POP 个个体
    void memoryCurrentPopulation() {

        nVec.clear();
        nVec = inVec;
        sort(nVec.begin(), nVec.end(), compByCost);
        for (unsigned int i = nVec.size() - 1; i > nVec.size() - N_POP - 1; i--) {
            nVec.pop_back();
        }

    }

    // 小生境遗传
    void nicheGA() {
        vector<Individual> mVec = inVec;
        int currentSize = N_POP + POP_SCALE;
        for (int i = 0; i < N_POP; i++) {
            mVec.push_back(nVec[i]);
        }
        for (int i = 0; i < currentSize - 1; i++) {
            for (int j = i + 1; j < currentSize; j++) {
                double difCount = sqrt(double((mVec[i].bitIn ^ mVec[j].bitIn).count()));
                if (difCount < L) {

                    if (mVec[i].cost < mVec[j].cost) {
                        mVec[j].fitness = 1e-30;
                    } else {
                        mVec[i].fitness = 1e-30;
                    }

                    //break;
                }
            }
        }
        sort(mVec.begin(), mVec.end(), compByFitness);
        inVec.clear();
        for (int i = 0; i < POP_SCALE; i++) {
            inVec.push_back(mVec[i]);
        }
    }

    // 产生新的一代
    void generateNextPopulation() {
        selectByLeague();
        //selectByRoulette();
        crossOverUniform();
        //crossOverOnePoint();
        mutationUniform();
        //mutationOnePoint();
    }

    // 整个进化的全过程
    void epoch() {

        gen = 1;
        // 生成初代
        generateInitalPopulation();
        // 计算 cost, fitness, 找出最优最差个体
        evalutePopulation();
        if (nodesNum > 200) { // 小case 不使用小生境
            memoryCurrentPopulation();
        }
        //show();
        while(gen < MAX_GENERATION) {
            gen++;
            // 通过选择 交叉 变异 生成下一代
            generateNextPopulation();
            mustChooseVec();
            // 评估新的一代
            evalutePopulation();
            if (nodesNum > 200) { // 小case 不使用小生境
                nicheGA();
                memoryCurrentPopulation();
            }

            // 精英策略
            //performEvolution();
            //show();

            double duration = double(clock() - start) / CLOCKS_PER_SEC;
            if (duration > ALL_TIME) {
                return;
            }
        }
    }

    // 输出种群现状
    void show() {
        int i;
        int sum = 0;
        double average;
        for (i = 0; i < POP_SCALE; i++) {
            sum += inVec[i].cost;
        }
        average = double((sum * 1.0) / (POP_SCALE * 1.0));
        printf("Gen = %d, Average Cost = %lf, Min Cost = %d, Current MinCost = %d, Current worstCost = %d\n",
               gen, average, everBestIndividual.cost, bestIndividual.cost, worstIndividual.cost);

/*
        vector<Individual>::iterator it;
        for (it = inVec.begin(); it != inVec.end(); it++) {
            (*it).show();
        }

        cout<<endl<<endl;
*/
    }
};
