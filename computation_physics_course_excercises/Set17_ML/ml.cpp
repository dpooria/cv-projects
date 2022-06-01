/*!
    Program for exercise set 17
    you can build and run this with main17.sh
*/

#include "data_info.h"
#include "svm.h"
#include <armadillo>
#include <bitset>
#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>
#include <sstream>

using namespace arma;

bool load_hoda(umat &x_train, umat &y_train, umat &x_test, umat &y_test);
bool load_boston(mat &xb_train, vec &yb_train, mat &xb_test, vec &yb_test);
umat bayes(const umat &x_train, const umat &y_train, const umat &x_test);
umat knn(const umat &x_train, const umat &y_train, const umat &x_test);
svm_model *m_svm_train(const umat &x_train, const umat &y_train);
umat svm_test(const svm_model *model, const umat &x_train);
mat confusionMat(const umat &y_obs, const umat &y_t);
template <int bit_size>
constexpr std::vector<std::bitset<bit_size>> getBitsetVector(const umat &dat);

int main(void)
{
    // Boston
    std::cout << "***Boston prices***" << std::endl;
    mat xb_tr, xb_ts;
    vec yb_tr, yb_ts;
    if (!load_boston(xb_tr, yb_tr, xb_ts, yb_ts))
    {
        std::cout << "error loading data" << std::endl;
        std::cout << "Please first run 'read_data.py'" << std::endl;
        return -1;
    }
    // fitting training data with Exact solution
    vec w_boston = inv(xb_tr.t() * xb_tr) * xb_tr.t() * yb_tr;
    vec y_lg = xb_ts * w_boston;
    double X2 = sum(abs(y_lg - yb_ts));
    std::ofstream fboston("boston-results.txt");
    for (uint i = 0; i < y_lg.size(); ++i)
    {
        fboston << yb_ts[i] << '\t' << y_lg[i] << '\n';
    }
    fboston.close();
    std::cout << "X^2 = " << X2 << std::endl;
    const char *gp_boston = R"E(gnuplot -e "
    set terminal png enhanced;set grid;
    set output 'boston-results.png'; set key bottom right;
    p 'boston-results.txt' u 1 w lp pt 6 ps 1.0 lc 'red' t 'actual',
        '' u 2 w lp pt 7 ps 1.0 lc 'blue' t 'predicted';")E";
    if (system(gp_boston))
    {
        std::cout << "plot error!" << std::endl;
    }
    // Hoda
    std::cout << "***Hoda data***" << std::endl;
    umat x_train, x_test, y_train, y_test;
    std::cout << "loading the data..." << std::endl;
    if (!load_hoda(x_train, y_train, x_test, y_test))
    {
        std::cout << "error loading data!" << std::endl;
        std::cout << "please first run 'read_data.py'" << std::endl;
        return -1;
    }
    std::clock_t ct = clock();
    const int n_test = x_test.n_rows;
    //===========================================================
    std::cout << "fitting with naive Bayes model..." << std::endl;
    umat y_bayes = bayes(x_train, y_train, x_test);
    std::cout
        << "Naive Bayes finished in "
        << double(std::clock() - ct) / CLOCKS_PER_SEC << " sec" << std::endl;
    auto C_bayes = confusionMat(y_bayes, y_test);
    // digonal elements of the confusion matrix are the fractions of right labeled features
    // so the sum of them times 100 is the percentege of right labeled ones
    std::cout << "Naive Bayes accuracy = "
              << sum(C_bayes.diag(0)) * 100 << "%" << std::endl;
    //===========================================================
    //svm and knn are slow and knn doesn't need that much of data so we just use 10% of the training data
    //for more accuracy increase the value of n_train_small
    uword n_train_small = std::lround(0.1 * x_train.n_rows);
    umat x_train_small = x_train.rows(0, n_train_small - 1);
    umat y_train_small = y_train.rows(0, n_train_small - 1);
    ct = std::clock();
    std::cout << "fitting with KNN method..." << std::endl;
    // umat y_knn = knn(x_train, y_train, x_test);
    umat y_knn = knn(x_train_small, y_train_small, x_test);
    std::cout << "KNN finished in "
              << double(std::clock() - ct) / CLOCKS_PER_SEC << " sec" << std::endl;
    auto C_knn = confusionMat(y_knn, y_test);
    std::cout << "KNN accuracy = "
              << sum(C_knn.diag(0)) * 100 << "%" << std::endl;
    //============================================================
    // SVM method is designed for two class classification problems
    // and in order to fit our multiclass problem libsvm uses a one vs one
    // algorithm. Therefore it will produce 45 ( 10 * 9 / 2) binary svm models(hence it is slow)
    ct = std::clock();
    std::cout << "***SVM***\nit'll take a while(about 4 minutes)" << std::endl;
    std::cout << "training the SVM model..." << std::endl;
    svm_model *model = m_svm_train(x_train_small, y_train_small);
    std::cout << "fitting with the SVM model..." << std::endl;
    umat y_svm = svm_test(model, x_test);
    svm_free_and_destroy_model(&model);
    std::cout << "SVM finished in "
              << double(std::clock() - ct) / CLOCKS_PER_SEC << " sec" << std::endl;
    auto C_svm = confusionMat(y_svm, y_test);
    std::cout << "SVM accuracy = " << sum(C_svm.diag(0)) * 100 << "%" << std::endl;
    //===============================================================
    // // save and plot the confusion matrices
    C_bayes.save("bayes-confusionMatrix.txt", arma::raw_ascii);
    C_knn.save("knn-confusionMatrix.txt", arma::raw_ascii);
    C_svm.save("svm-confusionMatrix.txt", arma::raw_ascii);
    const char *gnuplot = R"END(gnuplot -e "
    set terminal png enhanced;
    set xrange [-0.5:9.5]; set yrange [-0.5:9.5];
    set xtics 0,1,9; set ytics 0,1,9;
    set grid; set xlabel 'predicted'; set ylabel 'actual';
    set output 'bayes-confustionMatrix.png';
    p 'bayes-confusionMatrix.txt' matrix w image;
    set output 'knn-confustionMatrix.png';
    p 'knn-confusionMatrix.txt' matrix w image;
    set output 'svm-confustionMatrix.png';
    p 'svm-confusionMatrix.txt' matrix w image;"
    )END";
    if (system(gnuplot))
    {
        std::cout << "error plot!" << std::endl;
        std::cout << gnuplot << std::endl;
    }
    // plot 5 of the pictures randomly
    int pix = sqrt(m);
    for (int i = 1; i <= 5; ++i)
    {
        int a = randi(distr_param(0, n_test));
        umat pic = reshape(x_test.row(a), pix, pix);
        uvec actual = find(y_test.row(a) == 1);
        uvec knn = find(y_knn.row(a) == 1);
        uvec bayess = find(y_bayes.row(a) == 1);
        uvec svm = find(y_svm.row(a) == 1);
        // inplace_trans(f);
        pic.save("pic.txt", arma::raw_ascii);
        std::stringstream gp;
        gp << "gnuplot -e \"set terminal png;unset key;"
           << "set xrange [-1.25:32]; set yrange [-4.25:32];"
           << "unset border; unset xtics; unset ytics;"
           << "set output 'actual" << actual[0] << ".png';"
           << "set label 'bayes " << bayess[0] << "' at -1, -1.5;"
           << "set label 'knn " << knn[0] << "' at -1, -2.75;"
           << "set label 'svm " << svm[0] << "' at -1, -4;"
           << "p 'pic.txt' matrix w image;\"";
        if (system(gp.str().c_str()))
        {
            std::cout << "plot error!" << std::endl;
            std::cout << gp.str() << std::endl;
        }
    }
    return 0;
}

bool load_boston(mat &xb_train, vec &yb_train, mat &xb_test, vec &yb_test)
{
    arma::arma_rng::set_seed_random();
    mat x;
    vec y;
    if (!x.load(fnameBostonFeat) || !y.load(fnameBostonPrice))
    {
        return false;
    }
    const int nb = x.n_rows;
    uvec indx = linspace<uvec>(0, nb - 1, nb);
    indx = shuffle(indx);
    const int nt = std::lround(0.8 * nb);
    uvec ind_tr = indx.rows(0, nt - 1);
    uvec ind_ts = indx.rows(nt, nb - 1);
    xb_train = x.rows(ind_tr);
    yb_train = y(ind_tr);
    xb_test = x.rows(ind_ts);
    yb_test = y(ind_ts);
    return true;
}

/// Load data from files fnameX and fnameY (file names of armadillo binary files determined in data_info.h)
/// then split it into training and test datasets. you should have runned read_data.py(therefore read_data.cpp)
///  in order to this to work('.npz' data format should've been converted to arma_binary format)
/// @param[out] x_train -> train features from the file fnameX
/// @param[out] y_train -> train labels corresponding to x_train from the file fnameY
/// @param[out] x_test -> test features from the file fnameX
/// @param[out] y_test -> test labels corresponding to x_test from the file fnameY
/// @returns true if it succeeded in loading data or false when it fails to load
bool load_hoda(umat &x_train, umat &y_train, umat &x_test, umat &y_test)
{
    arma::arma_rng::set_seed_random();
    arma::Mat<int> x;
    arma::Mat<int> y;
    if (!x.load(fnameX) || !y.load(fnameY))
    {
        return false;
    }
    // split data into training(80%) and test(20%)
    const int n_train = std::lround(0.8 * n);
    uvec indx = linspace<uvec>(0, n - 1, n);
    // to choose training and test data randomly shuffle indexes and then
    //   chooe n_train and n_test elements from x
    indx = arma::shuffle(indx);
    uvec ind_train = indx.rows(0, n_train - 1);
    uvec ind_test = indx.rows(n_train, n - 1);
    x_train = conv_to<umat>::from(x.rows(ind_train));
    x_test = conv_to<umat>::from(x.rows(ind_test));
    y_train = conv_to<umat>::from(y.rows(ind_train));
    y_test = conv_to<umat>::from(y.rows(ind_test));
    return true;
}

/// Naiive Bayes model
/// @param x_train binary matrix of size (n_train x n_m) -> features with label
/// @param y_train binary matrix of size (n_train x n_c) -> labels corresponding to x_trian
/// @param x_test binary matrix of size (n_test x n_m) -> features without label
/// @returns binary matrix of size (n_test x n_c) -> predicted labels corresponding to x_test
umat bayes(const umat &x_train, const umat &y_train, const umat &x_test)
{
    const int n_train = x_train.n_rows;
    const int n_test = x_test.n_rows;
    const int n_c = y_train.n_cols;
    const int n_m = x_train.n_cols;
    // p_cond ~ p(x_i | c_j)
    mat p_cond(n_m, n_c, fill::zeros);
    // p_class ~ p(c)
    rowvec p_class = conv_to<rowvec>::from(sum(y_train, 0));
    p_class /= n_train;
    // count the non-zero elements in x_train and their label(class)
    for (int i = 0; i < n_train; ++i)
    {
        uvec cl = arma::find(y_train.row(i) == 1);
        uvec xi = arma::find(x_train.row(i) == 1);
        p_cond(xi, cl) = p_cond(xi, cl) + 1.0;
    }
    // normalize p_cond
    for (int i = 0; i < n_c; ++i)
    {
        p_cond.col(i) = (p_cond.col(i) + 1.0) / sum(p_cond.col(i) + 1.0);
    }
    p_cond = log(p_cond);
    umat y_res(n_test, n_c, fill::zeros);
    for (int i = 0; i < n_test; ++i)
    {
        // p the rows of p_cond corresponding to non-zeros elements of x_test(i)
        mat p = p_cond.rows(find(x_test.row(i) == 1));
        // p_c ~ p(Class | pic)
        rowvec p_c = prod(p, 0) % p_class;
        // the maximum index is the right label(class)
        y_res(i, p_c.index_max()) = 1;
    }
    return y_res;
}

/// KNN method(with 1 neighbor) and Hamming distance
/// @returns labels fitted to x_test
/// @param x_train a binary matrix of size (n_train x m) -> features with labels
/// @param y_train a binary matrix of size (n_train x c) -> corresponding labels of x_train
/// @param x_test a binary matrix of size (n_test x m) -> features without labels
umat knn(const umat &x_train, const umat &y_train, const umat &x_test)
{
    const int n_train = x_train.n_rows;
    const int n_test = x_test.n_rows;
    // convert x_train and x_test to std::vector<std::bitset<m>>
    // for efficiency purposes
    auto xtr = getBitsetVector<m>(x_train);
    auto xts = getBitsetVector<m>(x_test);
    arma::umat y_knn(n_test, c, fill::zeros);
    for (int i = 0; i < n_test; ++i)
    {
        int odmax = -999;
        int k = -1;
        for (int j = 0; j < n_train; ++j)
        {
            // compute the opposite Hamming distance(more is better)
            //  by counting the matching ones and matching zeros
            int od = (xts[i] & xtr[j]).count() + (~xts[i] & ~xtr[j]).count();
            // the index corresponding to maximum od (hence k ) is the closest point to x_test[i]
            if (od >= odmax)
            {
                odmax = od;
                k = j;
            }
        }
        y_knn.row(i) = y_train.row(k);
    }
    return y_knn;
}

/// function for knn to convert arma::umat to std::vector<std::bitset<bit_size>> class
/// @param bit_size integer size of bitset class -> constexpr
/// @param dat a binary matrix with <bit_size> cols
/// @returns convert of dat to std::vector<std::bitset<bit_size>>
template <int bit_size>
constexpr std::vector<std::bitset<bit_size>> getBitsetVector(const umat &dat)
{
    const int N = dat.n_rows;
    std::vector<std::bitset<bit_size>> dat_bit(N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < bit_size; ++j)
        {
            dat_bit[i][j] = dat(i, j);
        }
    }
    return dat_bit;
}

/// computes confusion matrix proportional to y_obs and y_t
mat confusionMat(const umat &y_obs, const umat &y_t)
{
    const int c_y = y_obs.n_cols;
    mat X(c_y, c_y, fill::zeros);
    for (uword i = 0; i < y_obs.n_rows; ++i)
    {
        // find the location of 1 in y_obs and y_t
        uvec r = find(y_obs.row(i) == 1);
        uvec c = find(y_t.row(i) == 1);
        // add one in the row and column proportional to the class
        // so if the both of them are the same,diagonal elements of matrix will increase
        X(r, c) = X(r, c) + 1.0;
    }
    X /= y_obs.n_rows;
    return X;
}

/// Interface for libsvm svm_train function
/// libsvm: https://www.csie.ntu.edu.tw/~cjlin/libsvm/
/// @returns svm_model struct obtained from svm_train function
svm_model *m_svm_train(const umat &x_train, const umat &y_train)
{
    // const int c_y = y_train.n_cols;
    const int n_train = x_train.n_rows;
    const int m_feat = x_train.n_cols;
    struct svm_parameter param;
    // svm_type -> multi class
    param.svm_type = C_SVC;
    // using polynomial kernel K(x, y) = (g*x'y + r)^d
    param.kernel_type = POLY;
    // degree -> d = 2 kerenel parameter
    param.degree = 2;
    // gamma -> g kernel parameter
    param.gamma = 1.0 / m_feat; // 1/num_features
    // coeff0 -> r kernel parameter
    param.coef0 = 0;
    param.cache_size = 100;
    param.C = 1; //the penalty parameter of the error term
    // other parameters are default parameters from svm-train.c
    param.eps = 1e-3;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 0;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;
    // prob -> svm_problem struct to save x and y and n_test
    struct svm_problem prob;
    // result model
    struct svm_model *model;
    prob.l = n_train;
    // y_labels -> converted y_train(labels) to libsvm format
    double *y_labels = new double[n_train];
    // x_spave -> converted x_train(features) to libsvm format
    struct svm_node **x_space = new svm_node *[n_train];
    for (int i = 0; i < n_train; ++i)
    {
        x_space[i] = new svm_node[m_feat + 1];
        for (int j = 0; j < m_feat; ++j)
        {
            // index -> feature label
            x_space[i][j].index = j;
            // value -> feature value
            x_space[i][j].value = x_train(i, j);
        }
        x_space[i][m_feat].index = -1; // to end the row
        // convert labels e.g. instead of {0,0,1} -> 2
        uvec cl = find(y_train.row(i) == 1);
        y_labels[i] = cl[0];
    }
    prob.x = x_space;
    prob.y = y_labels;
    // train the model
    model = svm_train(&prob, &param);
    // clean up
    for (int i = 0; i < n_train; ++i)
    {
        delete[] x_space[i];
    }
    delete[] x_space;
    delete[] y_labels;
    return model;
}

/// Interface function to libsvm svm_predict function
/// @param model svm_model struct obtained from libsvm svm_train finction
/// @returns predicted labels of x_test
umat svm_test(const svm_model *model, const umat &x_test)
{
    const int n_test = x_test.n_rows;
    const int m_feat = x_test.n_cols;
    umat y_svm(n_test, c, fill::zeros);
    for (int i = 0; i < n_test; ++i)
    {
        // find the on pixels
        uvec xi = find(x_test.row(i) != 0);
        uword nx = xi.size();
        // create a svm_node array of size of pixels that are on
        svm_node *x_space = new svm_node[nx + 1];
        for (int j = 0; j < nx; ++j)
        {
            // index is the pixel location
            x_space[j].index = xi[j];
            // value is 1 = on
            x_space[j].value = 1;
        }
        // end the array
        x_space[nx].index = -1;
        // make prediction with svm_predict function
        uword yc = (uword)svm_predict(model, x_space);
        // convert it back to our dataset format
        y_svm(i, yc) = 1;
        delete[] x_space;
    }
    return y_svm;
}
