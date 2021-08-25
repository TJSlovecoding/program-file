
#include <iostream>
#include<cmath>
#include<vector>
#include<time.h>

using namespace std;

int n = 0; //电梯数量
int m = 0; //楼层数量

double mpeople = 0.06;  //人的质量
double mcar = 0.8;  //电梯梯厢的质量
double g = 10;  //重力加速度
double EC = 0.0225;//每次启停的能耗

int limit = 0;//电梯限制人数

const int NC_max = 200;  //迭代次数
const double Alpha = 2;		//表征信息素重要程度的参数
const double Beta = 5;		//表征启发式因子重要程度的参数
const double Rho = 0.05;		//信息素蒸发系数
const double Q = 100000;		//信息素增加强度系数

double** h;//表示电梯之间的高度[m][m]

double** E;//表示从i楼到j楼消耗的能量[m][m]

double** Tau;//第i台电梯服务于第j层呼叫的信息素[n][m]
vector<int>* TABU;//存储电梯走过的路径[n]

vector<int>* updownrenshu;//保存每次上下的人数信息

int* uprenshu;//m层要上楼的人数[m]

int* downrenshu;//m层要下楼的人数[m]

int** renshu;//从a层到b层的人数[m][m]

int** liftrenshu;//电梯内还需要到达的人数[n][m]

vector<int>** bestTABU;//保存所有电梯路径信息[NC_max][n]
vector<int>** bestupdown;//保存每轮迭代的上下信息
double* bestE;//保留所有能耗信息[NC_max]

class people {
public:
    int start;
    int end;
    int dir;
};

void init() {

    cout << "请输入电梯个数" << endl;
    cin >> n;
    cout << "请输入楼层数：" << endl;
    cin >> m;
    cout << "请输入电梯限制人数：" << endl;
    cin >> limit;

    h = new double* [m];
    for (int i = 0; i < m; i++) {//二维动态数组
        h[i] = new double[m];
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (i != j) {
                h[i][j] = 3 * (double)abs(j - i);//每层楼三米
            }
            else {
                h[i][j] = 0;//设置不可达（原DBL_EPSILON）
            }
        }
    }

    E = new double* [m];
    for (int i = 0; i < m; i++) {//二维动态数组
        E[i] = new double[m];
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {//起初电梯里没有人
            E[i][j] = (mcar * g * h[i][j]) + 2 * EC;//一次启动一次停
        }
    }

    Tau = new double* [n];
    for (int i = 0; i < n; i++) {//二维动态数组
        Tau[i] = new double[m];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            Tau[i][j] = 1;
        }
    }

    uprenshu = new int[m];
    downrenshu = new int[m];
    for (int i = 0; i < m; i++) {
        uprenshu[i] = 0;
    }
    for (int j = 0; j < m; j++) {
        downrenshu[j] = 0;
    }

    renshu = new int* [m];
    for (int i = 0; i < m; i++) {//二维动态数组
        renshu[i] = new int[m];
    }

    for (int i = 0; i < m; i++) {//初始化人数
        for (int j = 0; j < m; j++) {
            renshu[i][j] = 0;
        }
    }

    liftrenshu = new int* [n];
    for (int i = 0; i < n; i++) {//二维动态数组
        liftrenshu[i] = new int[m];
    }

    for (int i = 0; i < n; i++) {//初始化每个电梯内的人数
        for (int j = 0; j < m; j++) {
            liftrenshu[i][j] = 0;
        }
    }

    TABU = new vector<int>[n];
    updownrenshu = new vector<int>[n];

    bestTABU = new vector<int>*[NC_max];
    for (int i = 0; i < NC_max; i++) {//二维动态数组
        bestTABU[i] = new vector<int>[n];
    }

    bestupdown = new vector<int>*[NC_max];
    for (int i = 0; i < NC_max; i++) {
        bestupdown[i] = new vector<int>[n];
    }

    bestE = new double[NC_max];
    for (int i = 0; i < NC_max; i++) {
        bestE[i] = 0.0;
    }
}

double rnd(double lower, double uper)	//生成lower和uper之间的一个double类型随机数
{
    return  (rand() / (double)RAND_MAX) * (uper - lower) + lower;
}

void ini2() {
    int numren = 0;
    cout << "请输入乘坐电梯的人的数量：" << endl;
    cin >> numren;
    vector<people>pp;
    srand((int)time(0));
    for (int i = 0; i < numren; i++) {//初始化乘电梯的人数
        people temp;
        temp.start = rand() % m;
        do {
            temp.end = rand() % m;
        } while (temp.end == temp.start);//避免在同一楼层
        if (temp.start < temp.end) {
            temp.dir = 1;//方向向上
        }
        else {
            temp.dir = 0;//方向向下
        }
        pp.push_back(temp);
    }

    for (int i = 0; i < pp.size(); i++) {
        renshu[pp[i].start][pp[i].end]++;
    }


    for (int i = 0; i < pp.size(); i++) {//初始化uprenshu和downrenshu
        if (pp[i].dir == 1) {
            uprenshu[pp[i].start]++;
        }
        else {
            downrenshu[pp[i].start]++;
        }
    }

}

void ant() {
    int NC = 0;

    int** renshu1 = new int* [m];//中间人数
    for (int i = 0; i < m; i++) {//二维动态数组
        renshu1[i] = new int[m];
    }

    int* uprenshu1 = new int[m];//中间人数

    int* downrenshu1 = new int[m];//中间人数
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            renshu1[i][j] = renshu[i][j];
        }
    }
    for (int j = 0; j < m; j++) {
        uprenshu1[j] = uprenshu[j];
    }
    for (int j = 0; j < m; j++) {
        downrenshu1[j] = downrenshu[j];
    }

    srand((int)time(0));

    int* a = new int[n];//存放各电梯的初始位置
    for (int i = 0; i < n; i++) {
        a[i] = rand() % m;//0~19
    }

    double* Ezong = new double[n];//这一代第n台电梯消耗的总能量
    int* currentrenshu = new int[n];

    while (NC < NC_max) {

        for (int i = 0; i < n; i++) {
            TABU[i].push_back(a[i]);//存储电梯初始位置
            updownrenshu[i].push_back(0);
        }

        for (int i = 0; i < n; i++) {
            Ezong[i] = 0;
        }

        for (int i = 0; i < n; i++) {
            currentrenshu[i] = 0;//当前电梯内的人数
        }

        int cc = 0;//判断变量
        do {
            cc = 0;
            for (int i = 0; i < m; i++) {//先向上服务

                if (uprenshu1[i] > 0) {//如果这一层有人要上

                    for (int a = 0; a < n; a++) {
                        for (int b = 0; b < m; b++) {
                            if (liftrenshu[a][b] > 0 && b <= i) {//如果这一层有人要下电梯
                                Ezong[a] += (E[TABU[a].back()][b] + currentrenshu[a] * mpeople * g * h[TABU[a].back()][b]);//移动到这一层的能量
                                currentrenshu[a] -= liftrenshu[a][b];
                                updownrenshu[a].push_back(liftrenshu[a][b]);//更新这一层出电梯人数
                                liftrenshu[a][b] = 0;//到这一层楼的人全出来了
                                TABU[a].push_back(b);//电梯移动到这一层
                            }
                        }
                    }

                    vector<double>p;
                    vector<int>xiabiao;//保存电梯的下标
                    double Psum = 0.0;		//概率值和
                    double rate = 0.0;		//随机数
                    double choose = 0.0;	//轮盘赌算法累加值
                    int yixuanlift = 0;//表示选中的电梯的编号
                    int choice = 0;//判断要不要停
                    for (int j = 0; j < n; j++) {//确定待选电梯
                        if (currentrenshu[j] < limit) {//如果还有剩位子
                            p.push_back(0.0);
                            xiabiao.push_back(j);
                            choice = 1;
                        }
                    }
                    if (choice == 1) {//只有有空闲电梯才服务这一层
                        for (int j = 0; j < p.size(); j++) {
                            double temp = 0;//中间能耗
                            temp = E[TABU[xiabiao[j]].back()][i] + currentrenshu[xiabiao[j]] * mpeople * g * h[TABU[xiabiao[j]].back()][i];
                            p[j] = pow(Tau[xiabiao[j]][i], Alpha) * pow((1 / temp), Beta);
                            Psum += p[j];//概率和
                        }

                        rate = rnd(0.0, Psum);
                        for (int k = 0; k < p.size(); k++)//赌轮盘算法选电梯
                        {
                            choose += p[k];
                            if (choose > rate)//选到了
                            {
                                yixuanlift = xiabiao[k];
                                break;
                            }
                        }

                        Ezong[yixuanlift] += (E[TABU[yixuanlift].back()][i] + currentrenshu[yixuanlift] * mpeople * g * h[TABU[yixuanlift].back()][i]);//更新总消耗能量
                        int temp1 = currentrenshu[yixuanlift];//保存有人上电梯前的电梯人数
                        for (int a = i + 1; a < m; a++) {//更新电梯内的人数和等候人数
                            int temp = limit - currentrenshu[yixuanlift];//剩余席位
                            if (temp > 0) {
                                if (temp >= renshu1[i][a]) {//如果电梯剩余人数够这一层的人进来
                                    liftrenshu[yixuanlift][a] += renshu1[i][a];//这一层所有要向上的人都进电梯
                                    uprenshu1[i] -= renshu1[i][a];//更新这一层要上楼的人
                                    currentrenshu[yixuanlift] += renshu1[i][a];
                                    renshu1[i][a] = 0;

                                }
                                else {
                                    liftrenshu[yixuanlift][a] += temp;
                                    uprenshu1[i] -= temp;
                                    renshu1[i][a] -= temp;
                                    currentrenshu[yixuanlift] += temp;
                                }
                            }
                            else {
                                break;
                            }
                        }
                        int temp2 = currentrenshu[yixuanlift] - temp1;
                        updownrenshu[yixuanlift].push_back(temp2);//更新这一层上电梯的人数
                        TABU[yixuanlift].push_back(i);//保存已经走的路径
                    }

                }
            }

            for (int a = 0; a < n; a++) {//把电梯内剩下的人都送出去(清零)
                for (int b = 0; b < m; b++) {
                    if (liftrenshu[a][b] > 0) {//如果这一层有人要下电梯

                        Ezong[a] += (E[TABU[a].back()][b] + currentrenshu[a] * mpeople * g * h[TABU[a].back()][b]);//移动到这一层的能量
                        currentrenshu[a] -= liftrenshu[a][b];
                        updownrenshu[a].push_back(liftrenshu[a][b]);//更新出电梯的人数

                        liftrenshu[a][b] = 0;//到这一层楼的人全出来了

                        TABU[a].push_back(b);//电梯移动到这一层
                    }
                }
            }
            //分界线hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
            for (int i = m - 1; i >= 0; i--) {//再向下服务

                if (downrenshu1[i] > 0) {//如果这一层有人要下

                    for (int a = 0; a < n; a++) {
                        for (int b = m - 1; b >= 0; b--) {
                            if (liftrenshu[a][b] > 0 && b >= i) {//如果这一层有人要下电梯
                                Ezong[a] += (E[TABU[a].back()][b] + currentrenshu[a] * mpeople * g * h[TABU[a].back()][b]);//移动到这一层的能量
                                currentrenshu[a] -= liftrenshu[a][b];
                                updownrenshu[a].push_back(liftrenshu[a][b]);//更新出电梯的人数
                                liftrenshu[a][b] = 0;//到这一层楼的人全出来了
                                TABU[a].push_back(b);//电梯移动到这一层
                            }
                        }
                    }

                    vector<double>p;
                    vector<int>xiabiao;//保存电梯的下标
                    double Psum = 0.0;		//概率值和
                    double rate = 0.0;		//随机数
                    double choose = 0.0;	//轮盘赌算法累加值
                    int yixuanlift = 0;//表示选中的电梯的编号
                    int choice = 0;//判断要不要停
                    for (int j = 0; j < n; j++) {//确定待选电梯
                        if (currentrenshu[j] < limit) {
                            p.push_back(0.0);
                            xiabiao.push_back(j);
                            choice = 1;
                        }
                    }
                    if (choice == 1) {//只有有空闲电梯才服务这一层
                        for (int j = 0; j < p.size(); j++) {
                            double temp = 0;//中间能耗
                            temp = E[TABU[xiabiao[j]].back()][i] + currentrenshu[xiabiao[j]] * mpeople * g * h[TABU[xiabiao[j]].back()][i];
                            p[j] = pow(Tau[xiabiao[j]][i], Alpha) * pow((1 / temp), Beta);
                            Psum += p[j];//概率和
                        }

                        rate = rnd(0.0, Psum);
                        for (int k = 0; k < p.size(); k++)//赌轮盘算法选电梯
                        {
                            choose += p[k];
                            if (choose > rate)//选到了
                            {
                                yixuanlift = xiabiao[k];
                                break;
                            }
                        }

                        Ezong[yixuanlift] += (E[TABU[yixuanlift].back()][i] + currentrenshu[yixuanlift] * mpeople * g * h[TABU[yixuanlift].back()][i]);//更新总消耗能量
                        int temp1 = currentrenshu[yixuanlift];
                        for (int a = i - 1; a >= 0; a--) {//更新电梯内的人数和等候人数
                            int temp = limit - currentrenshu[yixuanlift];//剩余席位
                            if (temp > 0) {
                                if (temp >= renshu1[i][a]) {//如果电梯剩余人数够这一层的人进来
                                    liftrenshu[yixuanlift][a] += renshu1[i][a];//这一层所有要向上的人都进电梯
                                    downrenshu1[i] -= renshu1[i][a];//更新这一层要上楼的人
                                    currentrenshu[yixuanlift] += renshu1[i][a];
                                    renshu1[i][a] = 0;

                                }
                                else {
                                    liftrenshu[yixuanlift][a] += temp;
                                    downrenshu1[i] -= temp;
                                    renshu1[i][a] -= temp;
                                    currentrenshu[yixuanlift] += temp;
                                }
                            }
                            else {
                                break;
                            }
                        }
                        int temp2 = currentrenshu[yixuanlift] - temp1;
                        updownrenshu[yixuanlift].push_back(temp2);
                        TABU[yixuanlift].push_back(i);//保存已经走的路径
                    }

                }
            }

            for (int a = 0; a < n; a++) {//把电梯内剩下的人都送出去
                for (int b = m - 1; b >= 0; b--) {
                    if (liftrenshu[a][b] > 0) {//如果这一层有人要下电梯
                        Ezong[a] += (E[TABU[a].back()][b] + currentrenshu[a] * mpeople * g * h[TABU[a].back()][b]);//移动到这一层的能量
                        currentrenshu[a] -= liftrenshu[a][b];
                        updownrenshu[a].push_back(liftrenshu[a][b]);
                        liftrenshu[a][b] = 0;//到这一层楼的人全出来了
                        TABU[a].push_back(b);//电梯移动到这一层
                    }
                }
            }

            for (int i = 0; i < m; i++) {//如果请求未结束
                if (uprenshu1[i] > 0) {//还有人要上楼
                    cc = 1;
                }
            }
            for (int i = 0; i < m; i++) {//如果请求未结束
                if (downrenshu1[i] > 0) {//还有人要下楼
                    cc = 1;
                }
            }

        } while (cc == 1);

        //更新信息素👇
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                Tau[i][j] = (1 - Rho) * Tau[i][j];//挥发系数
            }
        }

        for (int i = 0; i < n; i++) {
            for (int j = 1; j < TABU[i].size(); j++) {//初始层不属于服务层
                if (TABU[i][j] != TABU[i][j - 1]) {//防止同一层楼多次更新
                    Tau[i][j] += Q / Ezong[i];
                }
            }
        }

        //保存每次的总耗能
        for (int i = 0; i < n; i++) {
            bestE[NC] += Ezong[i];
        }
        //保存每次的路径
        for (int i = 0; i < n; i++) {
            bestTABU[NC][i] = TABU[i];
        }

        for (int i = 0; i < n; i++) {//保存每部电梯的进出人数数据
            bestupdown[NC][i] = updownrenshu[i];
        }

        for (int i = 0; i < n; i++) {//清空禁忌表
            TABU[i].clear();
        }
        for (int i = 0; i < n; i++) {//清空禁忌表
            updownrenshu[i].clear();
        }

        cout << bestE[NC] << endl;//查看当代的能量总和

        NC++;

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                renshu1[i][j] = renshu[i][j];
            }
        }
        for (int j = 0; j < m; j++) {
            uprenshu1[j] = uprenshu[j];
        }
        for (int j = 0; j < m; j++) {
            downrenshu1[j] = downrenshu[j];
        }
    }

    //打印出耗能最少的方案
    double minE = bestE[0];//保存最小能耗
    int minxiabiao = 0;//保存最小能耗的下标

    for (int i = 0; i < NC_max; i++) {
        if (bestE[i] < minE) {
            minE = bestE[i];
            minxiabiao = i;
        }
    }

    //cout << "最小的能耗为：" << bestE[minxiabiao] << endl;
    
    cout << "最小的能耗为：" << bestE[minxiabiao] << endl;
    for (int i = 0; i < n; i++) {
        cout << "电梯" << i + 1 << "的路径为: ";
        for (int j = 0; j < bestTABU[minxiabiao][i].size(); j++) {
            cout << bestTABU[minxiabiao][i][j] + 1 << " ";
        }
        cout << endl;
        cout << "电梯" << i + 1 << "的进出人数数据为：";
        for (int j = 0; j < bestupdown[minxiabiao][i].size(); j++) {
            cout << bestupdown[minxiabiao][i][j] << " ";
        }
        cout << endl;
    }

    //释放动态数组
    delete[]uprenshu1;
    delete[]downrenshu1;
    delete[]a;
    delete[]Ezong;
    delete[]currentrenshu;
    for (int i = 0; i < m; i++) {
        delete[]renshu1[i];
    }
    delete[]renshu1;

}

void Delete1() {//释放内存
    for (int i = 0; i < m; i++) {
        delete[]h[i];
    }
    delete[]h;

    delete[]E;

    delete[]Tau;

    delete[]TABU;
    delete[]uprenshu;
    delete[]downrenshu;
    for (int i = 0; i < m; i++) {
        delete[]renshu[i];
    }
    delete[]renshu;

    for (int i = 0; i < n; i++) {
        delete[]liftrenshu[i];
    }
    delete[]liftrenshu;

    delete[]bestTABU;
    delete[]bestupdown;
    delete[]updownrenshu;
    delete[]bestE;
}

int main() {
    init();
    ini2();
    ant();
    Delete1();
    return 0;

}