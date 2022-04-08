#include <iostream>
#include <random>
#include <math.h>

long long a;
long long c;
long long m;
long long seed;
long long prev_xn;

void init_rand(long long rand_seed, long long rand_a, long long rand_c, long long rand_m){
    seed = rand_seed;
    prev_xn = seed;
    a = rand_a;
    c = rand_c;
    m = rand_m;
}

double custom_rand(){
    prev_xn = (a * prev_xn + c) % m;
    return (double)prev_xn / (m-1);
}

double std_rand(){
    return (double)rand()/RAND_MAX;
}

void allocMatrix(int n, double** &mat){
    mat = (double**)malloc(sizeof(double*) * n);
    for(int i = 0; i < n; i++){
        mat[i] = (double*)malloc(sizeof(double) * 2);
    }
}

void freeMatrix(int n, double** &mat){
    for(int i = 0; i < n; i++){
        std::free(mat[i]);
    }
    std::free(mat);
}

int* calculatePi(int n, double** pointsIn, double** pointsOutside, double (*randFunc)()){
    int inCount = 0;
    int outCount = 0;
    for(int i = 0; i < n; i++) {
        double x = randFunc();
        double y = randFunc();
        if(sqrt(x*x + y*y) < 1){
            pointsIn[inCount][0] = x;
            pointsIn[inCount][1] = y;
            inCount++;
        } else{
            pointsOutside[outCount][0] = x;
            pointsOutside[outCount][1] = y;
            outCount++;
        }
    }
    int* out = (int*)malloc(sizeof(int)*2);
    out[0] = inCount;
    out[1] = outCount;
    return out;
}

void appendToFile(FILE* &fp, int n, double** mat){
    for(int i = 0; i < n; i++) {
        fprintf(fp, "%f %f\n", mat[i][0], mat[i][1]);
    }
    fprintf(fp, "\n\n");
}

void getGif(std::string fileName){
    FILE* gnupipe = _popen("gnuplot -persistent", "w");
    fprintf(gnupipe, "f(x)=sqrt(1-x*x)\n");
    fprintf(gnupipe, "cd \"C:\\\\Users\\\\Dinislam\\\\CLionProjects\\\\HPC_HW3\\\\cmake-build-debug\"\n");
    fprintf(gnupipe, ("stats \"" + fileName + "\" name \"A\" nooutput \n").c_str());
    fprintf(gnupipe, "set xrange[0:1]\n");
    fprintf(gnupipe, "set yrange[0:1]\n"
                            "set size ratio -1\n");
    fprintf(gnupipe, "set term gif animate delay 50 size 800, 800 crop\n  ");
    fprintf(gnupipe, ("set output \"" + fileName + ".gif\"\n").c_str());
    fprintf(gnupipe, ("do for [i=0:int((A_blocks-1)/2)] \\\n"
                     "{ stats \"" + fileName + "\" index i*2 name \"B\" nooutput \n"
                     "  stats \"" + fileName + "\" index i*2+1 name \"C\" nooutput \n"
                     "set title sprintf(\"n=%%.0f, Pi=%%.4f\", (i+1)*1000, B_records*4./(B_records+C_records)) \n" //" \"n=\".((i+1)*1000).\", Pi=\".sprintf(\"%f\", B_records*4./(B_records+C_records)) \n"
                     "plot \"" + fileName + "\" index i*2 u 1:2 pt 7 ps 0.15 notitle, \\\n"
                     "\""+fileName+"\" index i*2+1 u 1:2 pt 7 ps 0.15 notitle, \\\n"
                     "f(x) w lines lw 1.5 notitle}\n").c_str());
    fprintf(gnupipe, "exit");
    _pclose(gnupipe);
}

void getAccuracyGraph(const std::string& fileName1, const std::string& fileName2){
    FILE* gnupipe = _popen("gnuplot -persistent", "w");
    fprintf(gnupipe, "cd \"C:\\\\Users\\\\Dinislam\\\\CLionProjects\\\\HPC_HW3\\\\cmake-build-debug\"\n");
    fprintf(gnupipe, "set title \"PI estimation accuracy\"\n");
    fprintf(gnupipe, "f(x) = pi\n");
    fprintf(gnupipe, ("plot \"" + fileName1 + ".acc\" u 1:2 w lines title \"Built-in rand acc\", \\\n"
                      "\"" + fileName2 + ".acc\" u 1:2 w lines title \"Custom rand acc\", \\\n"
                      "f(x) w lines lw 2 title \"PI\"\n").c_str());
    fprintf(gnupipe, "exit");
    _pclose(gnupipe);
}

void getPiAccuracy(int start, int end, int step, const std::string& fileName, double (*randFunc)()){
    FILE* fp1;
    fp1 = fopen(fileName.c_str(), "w");
    FILE* fp2;
    fp2 = fopen((fileName + ".acc").c_str(), "w");

    for(int nPoints = start; nPoints <= end; nPoints += step){
        double** in;
        double** out;
        allocMatrix(nPoints, in);
        allocMatrix(nPoints, out);
        int* vals = (int*)malloc(sizeof(int)*2);
        vals = calculatePi(nPoints, in, out, randFunc);
        appendToFile(fp1, vals[0], in);
        appendToFile(fp1, vals[1], out);
        fprintf(fp2, "%d %f\n", nPoints, (double)vals[0]*4/nPoints);
        freeMatrix(vals[0], in);
        freeMatrix(vals[1], out);
        std::free(vals);
    }
}

double f(double x){
    return abs(sin(x)) + 5*exp(-pow(x,4))*cos(x);
}

double monteCarloIntegral(int n, double a, double b, double (*randFunc)(), double (*integratedFunc)(double)){
    double N = 1/(double)n;
    double interval = b-a;
    double ans = 0;
    for(int i = 0; i < n; i++){
        double randX = a + interval * randFunc();
        ans += integratedFunc(randX)*interval*N;
    }
    return ans;
}

void getIntegralAccuracy(
        int start, int end, int step,
        double a, double b,
        double (*randFunc)(), double (*integratedFunc)(double),
        std::string &fileName
        )
{
    FILE* fp = fopen(fileName.c_str(), "w");
    double ans;
    for(int nIter = start; nIter <= end; nIter *= step){
        ans = monteCarloIntegral(nIter, a, b, randFunc, integratedFunc);
        fprintf(fp, "%d %f\n", nIter, ans);
    }

}

void getIntegralGraph(std::string &fileName1, std::string &fileName2){
    double ans = 7.09544;
    FILE* gnupipe = _popen("gnuplot -persistent", "w");
    fprintf(gnupipe, "f(x) = 7.09544\n");
    fprintf(gnupipe, "cd \"C:\\\\Users\\\\Dinislam\\\\CLionProjects\\\\HPC_HW3\\\\cmake-build-debug\"\n");
    fprintf(gnupipe, "set title \"Accuracy of integration techniques\"\n");
    fprintf(gnupipe, "set xrange[0:1048576]\n");
    fprintf(gnupipe, "set logscale x 2\n");
    fprintf(gnupipe, "set xlabel \"Iterations(log scale)\"\n");
    fprintf(gnupipe, "set ylabel \"Result\"\n");
    fprintf(gnupipe, ("plot \"" + fileName1 + "\" u 1:2 w lines title \"Built-in random\", \\\n"
                      "\"" + fileName2 + "\" u 1:2 w lines title \"Custom random\", \\\n"
                      "f(x) w lines title \"Answer\"\n").c_str());
    fprintf(gnupipe, "exit");
    _pclose(gnupipe);
}

int main() {
    init_rand(1, 845, 2625, 8192);
    //std::string f1_name = "task_1_std_rand.tmp";
    //getPiAccuracy(1000, 30000, 1000, f1_name, std_rand);
    //getGif(f1_name);
    //std::string f2_name = "task_1_custom_rand.tmp";
    //getPiAccuracy(1000, 30000, 1000, f2_name, custom_rand);
    //getGif(f2_name);
    //getAccuracyGraph(f1_name, f2_name);

    //double stdRandAns = monteCarloIntegral(100000, 0, 5, std_rand, f);
    //double customRandAns = monteCarloIntegral(100000, 0, 5, custom_rand, f);

    std::string f3_name = "task_2_std_rand.tmp";
    getIntegralAccuracy(1024, (int)pow(2, 20), 2, 0, 5, std_rand, f, f3_name);

    std::string f4_name = "task_2_custom_rand.tmp";
    getIntegralAccuracy(1024, (int)pow(2, 20), 2, 0, 5, custom_rand, f, f4_name);

    getIntegralGraph(f3_name, f4_name);
    return 0;
}
