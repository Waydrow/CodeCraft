#include "deploy.h"
#include "lib/lib_time.h"
#include <cmath>

#define MAX_NODES_NUM 1300 // 最大网络节点数量

int MAX_GENERATION; // 迭代次数

double ALL_TIME = 87; // 程序最长运行时间限制 单位 s

int POP_SCALE; // 种群规模

int N_POP = 100; // 记忆前 N 个当前种群个体

double L = 3.2; // 最小 Hamming 距离

double Pc = 0.7; // 交叉概率
// 0.0012
double Pm = 0.0012; // 变异概率

double Ps = 0.8;

const double A = 9.903438;
double Pcmax = 0.8;
double Pcmin = 0.7;
double Pmmax = 0.0014;
double Pmmin = 0.0007;

double c = 1.75; // 线性变换参数
double m; // 指数变换参数

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
    vector<int> directVec;
    vector<int> notChooseVec;
    MinCostFlowSolution *mcf;

    ~Population() {}
public:
    Population(MinCostFlowSolution *mcfs) { // constructor
        mcf = mcfs;
        if (nodesNum < 800) {
            POP_SCALE = 300;
            MAX_GENERATION = 1000;
            Pm = 0.001;
            N_POP = 125;
            L = 1.4;
        } else {
            POP_SCALE = 78;
            MAX_GENERATION = 1000;
            Pm = 0.0008;
            N_POP = 19;
            L = 1.8;
        }
        // MAX_GENERATION = 2;
        // 设置指数变换参数
        m = 1 + log10(MAX_GENERATION);
        for (int i = 0; i < POP_SCALE; i++) {
            Individual ind = Individual(); // 初始化
            inVec.push_back(ind);
        }

        directVec = mcf->getMustChoose();
        //notChooseVec = mcf->getNotChoose();
        //sort(notChooseVec.begin(), notChooseVec.end());
    }

    // 随机选取第一代
    void generateInitalPopulation() {
        int i = 0;
        while (i < POP_SCALE / 2) { // 每个体
            // 随机产生服务器的个数
            int directSize = directVec.size();
            int serverNum = get_random_int(1, clientNum - directSize);
            set<int> s;
            for (int j = 0; j < directSize; j++) {
                s.insert(directVec[j]);
            }
            /*
            vector<int> nodesTmp;
            vector<int>::iterator itt = notChooseVec.begin();
            for (int j = 0; j < nodesNum; j++) {
                if (j == (*itt)) {
                    itt++;
                    continue;
                }
                nodesTmp.push_back(j);
            }
            */
            // 产生 serverNum 个随机数
            while (true) {
                int r = get_random_int(0, nodesNum - 1);
                s.insert(r);
                if (int(s.size()) == serverNum + directSize) {
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

    void generateInitialBetter() {
        for (int i = POP_SCALE / 2; i < POP_SCALE; i++) {
            for (int j = 0; j < nodesNum; j++) {
                double pp = get_random_real(0, 1);
                if (pp < nodesP[j])
                    inVec[i].bitIn.set(j);
            }
        }
    }
    // 计算 cost, 与适应度指数变换对应
    void calcCostValueForExp() {
        double sumFit = 0;
        for (int i = 0; i < POP_SCALE; i++) {

            pair<pair<long long,long long>, bitset<BITSIZE> > pair1 = mcf->CalCost(inVec[i].bitIn, 0, true);
            inVec[i].cost = pair1.first.first;

            sumFit += double(inVec[i].cost);

            if (inVec[i].cost < everBestIndividual.cost) {
                everBestIndividual = inVec[i];
            }

            if ((nodesNum < 800 && gen <= 35) || (nodesNum > 1000 && gen <= 65)) {
                inVec[i].bitIn = pair1.second;
                inVec[i].cost = pair1.first.second;
            }
        }
        averageFit = sumFit / (POP_SCALE * 1.0);
    }

    // 计算适应度, 指数变换
    void calcFitnessValueExp() {
        // 指数
        double aa = pow(gen, 1 / m) / (averageFit + 0.001);
        maxFit = 0;
        minFit = 1e10;
        double sumF = 0;
        for (int i = 0; i < POP_SCALE; i++) {
            inVec[i].fitness = exp(-1 * aa * double(inVec[i].cost));
            sumF += inVec[i].fitness;
            if (inVec[i].fitness < minFit) {
                minFit = inVec[i].fitness;
            }
            if (inVec[i].fitness > maxFit) {
                maxFit = inVec[i].fitness;
            }
        }
        averageFit = sumF / (POP_SCALE * 1.0);
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
    void directVecChoose() {
        int num = directVec.size();
        for (int i = 0; i < POP_SCALE; i++) {
            for (int j = 0; j < num; j++) {
                inVec[i].bitIn.set(directVec[j]);
            }
        }
    }

    // not choose
    void notChoose() {
        int num = notChooseVec.size();
        for (int i = 0; i < POP_SCALE; i++) {
            for (int j = 0; j < num; j++) {
                inVec[i].bitIn.reset(notChooseVec[j]);
            }
        }
    }

    // evaluate
    void evalutePopulation() {
        calcCostValueForExp();
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

    void crossOverTwoPoint() {
        int i, j;
        vector<int> index;
        int point, point2, temp;
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
                point = get_random_int(1, nodesNum-2);
                point2 = get_random_int(1, nodesNum-1);
                int lower = min(point, point2);
                int upper = max(point, point2);
                bita = inVec[index[i]].bitIn;
                bitb = inVec[index[i+1]].bitIn;
                for (j = 0; j < lower; j++) {
                    temp = bita[j];
                    bita[j] = bitb[j];
                    bitb[j] = temp;
                }
                for (j = upper; j < nodesNum; j++) {
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
/*
            double biggerFit = max(inVec[index[i]].fitness, inVec[index[i+1]].fitness);
            if (biggerFit >= averageFit) {
                Pc = Pcmin + (Pcmax - Pcmin) / (1 + exp(A * (2 * (biggerFit - averageFit) / (maxFit - averageFit))));
            } else {
                Pc = Pcmax;
            }
*/
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
                point = get_random_int(1, nodesNum-2);
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

    void crossOverFromBetter() {
        bitset<BITSIZE> bita, bitb;
        for (int i = 0; i < POP_SCALE; i++) {
            int pos = 0;
            int pos1 = get_random_int(0, POP_SCALE-1);
            int pos2 = get_random_int(0, POP_SCALE-1);
            if (inVec[pos1].cost < inVec[pos2].cost) {
                pos = pos1;
            } else {
                pos = pos2;
            }

            int poss = 0;
            int pos3 = get_random_int(0, POP_SCALE-1);
            int pos4 = get_random_int(0, POP_SCALE-1);
            if (inVec[pos3].cost < inVec[pos4].cost) {
                poss = pos3;
            } else {
                poss = pos4;
            }

            double p = get_random_real(0, 1);
            if (p < Pc) {
                bita = inVec[pos].bitIn;
                bitb = inVec[pos4].bitIn;
                for (int j = 0; j < nodesNum; j++) {
                    int t = get_random_int(0, 1);
                    if (t == 1) {
                        int temp = bita[j];
                        bita[j] = bitb[j];
                        bitb[j] = temp;
                    }
                }
                // 判断服务器数量是否超过消费节点数量, 若超过, 则不交叉
                if (bita.count() > 0 && bita.count() <= unsigned(clientNum) && bitb.count() > 0 && bitb.count() <= unsigned(clientNum)) {
                    inVec[pos].bitIn = bita;
                    inVec[poss].bitIn = bitb;
                }
            }
        }
    }

    // 均匀变异
    void mutationUniform() {
        double p;
        // bit mutation
        for (int i = 0; i < POP_SCALE; i++) {
/*
            if (inVec[i].fitness >= averageFit) {
                Pm = Pmmin + (Pmmax - Pmmin) / (1 + exp(A * (2 * (inVec[i].fitness - averageFit) / (maxFit - averageFit))));
            } else {
                Pm = Pmmax;
            }
*/
            for (int j = 0; j < nodesNum; j++) {
                p = get_random_real(0, 1);
                // mutation
                if (p < Pm) {
                    // 该位取反
                    inVec[i].bitIn.flip(j);
                    // 变异后服务器节点数量大于消费节点数量
                    if (inVec[i].bitIn.count() > unsigned(clientNum) || inVec[i].bitIn.count() == 0) {
                        inVec[i].bitIn.flip(j);
                    }
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
                if (inVec[i].bitIn.count() > unsigned(clientNum) && inVec[i].bitIn.count() > 0) {
                    inVec[i].bitIn.flip(postion);
                }
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
                        mVec[j].fitness = 1e-10;
                    } else {
                        mVec[i].fitness = 1e-10;
                    }
                }
            }
        }
        sort(mVec.begin(), mVec.end(), compByFitness);
        inVec.clear();
        for (int i = 0; i < POP_SCALE; i++) {
            inVec.push_back(mVec[i]);
        }
    }

    void swapGen() {
        for (int i = 0; i < POP_SCALE; i++) {
            double pp = get_random_real(0, 1);
            if (pp < Ps) {
                int n = 1;
                while(n--) {
                    int pos1 = get_random_int(0, nodesNum-1);
                    //int pos2 = get_random_int(0, nodesNum-1);
                    int pos2;
                    if (pos1 == nodesNum-1)
                        pos2 = pos1 - 1;
                    else
                        pos2 = pos1 + 1;
                    int t = inVec[i].bitIn[pos1];
                    inVec[i].bitIn[pos1] = inVec[i].bitIn[pos2];
                    inVec[i].bitIn[pos2] = t;
                }
            }
        }
    }

    // 产生新的一代
    void generateNextPopulation() {
        selectByLeague();
        // selectByRoulette();
        // crossOverTwoPoint();
        crossOverUniform();
        // crossOverOnePoint();
        // crossOverFromBetter();
        mutationUniform();
        // mutationOnePoint();
        // swapGen();
    }

    // 整个进化的全过程
    void epoch() {
        gen = 1;
        // 生成初代
        generateInitalPopulation();
        generateInitialBetter();
        // 计算 cost, fitness, 找出最优最差个体
        evalutePopulation();
        memoryCurrentPopulation();
        // show();
        while(gen < MAX_GENERATION) {
            gen++;
            // 通过选择 交叉 变异 生成下一代
            generateNextPopulation();
            // must choose direct node
            if (nodesNum > 1000) {
                directVecChoose();
                //notChoose();
            }
            // 评估新的一代
            evalutePopulation();
            nicheGA();
            memoryCurrentPopulation();
            // 精英策略
            performEvolution();
            // show();
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
        i = 0;
        for (it = inVec.begin(); it != inVec.end(); it++) {
            printf("%d:  ", i++);
            (*it).show();
        }

        cout<<endl<<endl;
*/
    }
};
