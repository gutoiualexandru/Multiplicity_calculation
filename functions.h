#include <iostream>
#include <vector>
#include <cmath>
#include <utility>
#include <limits>
#include <ROOT/TThreadExecutor.hxx>

#include "generators.h"
using namespace std;
/// ////////////////////////////////////////////////////////////
void draw_distribution(vector<double> distribution)
{
    vector<double> X;
    for (double i = 0; i < (double)distribution.size(); i++)
    {
        X.push_back(i);
    }
    auto *t1 = new TGraph(X.size(), X.data(), distribution.data());
    t1->Draw();
}
void print_distribution(vector<double> distribution)
{
    for (int i = 0; i < distribution.size(); i++)
    {
        cout << distribution[i];
        if (i != distribution.size() - 1)
        {
            cout << ",";
        }
    }
    cout << endl;
}
void product_vector_scalar(vector<double> &grad, double lr)
{
    for (double &element : grad)
    {
        element *= lr;
    }
}
void sum_vectors(vector<double> &grad, vector<double> &angles)
{
    for (size_t i = 0; i < angles.size(); ++i)
    {
        angles[i] += grad[i];
    }
}
// Function to compute factorial
void normalize(std::vector<double> &v)
{
    double sum = 0;
    for (double i : v)
    {
        sum += i;
    }
    for (double &i : v)
    {
        i /= sum;
    }
}
double factorial(int n)
{
    if (n == 0 || n == 1)
    {
        return 1;
    }
    double result = 1;
    for (int i = 2; i <= n; ++i)
    {
        result *= i;
    }
    return result;
}

// Function to compute the normalization constant c_m
double compute_c_m(int m, int N, double epsilon)
{
    double c_m = 0.0;
    for (int i = 0; i <= m; ++i)
    {
        double numerator = factorial(m) * factorial(N);
        double denominator = factorial(i) * factorial(m - i) * factorial(N - i);
        double term = numerator / denominator * std::pow(epsilon / N, i) * std::pow(1 - epsilon, m - i);
        c_m += term;
    }
    return c_m;
}

// Function to compute f(n, m, epsilon)
double f(int n, int m, int N, double epsilon)
{
    if (n > m)
    {
        return 0;
    }
    double c_m = compute_c_m(m, N, epsilon);

    double numerator = factorial(m) * factorial(N);
    double denominator = factorial(n) * factorial(m - n) * factorial(N - n);
    double term = numerator / denominator * std::pow(epsilon / N, n) * std::pow(1 - epsilon, m - n);

    return (1.0 / c_m) * term;
}
class Setup
{
public:
    // Function to set limits and initialize the matrix with m x n elements
    void set_limits(int m, int n)
    {
        this->m = m;
        this->n = n;
        matrix.resize(m, std::vector<double>(n, 0.0)); // Initialize all elements to 0.0
    }

    // Function to set epsilon0
    void set_epsilon0(double epsilon)
    {
        epsilon0 = epsilon;
    }

    // Function to set N
    void set_N(int N)
    {
        this->N = N;
    }

    // Function to fill the matrix with the values of f(n, m, epsilon0)
    void fill_matrix()
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                matrix[i][j] = f(i + 1, j + 1, N, epsilon0);
            }
        }
    }
    // Function to display the matrix
    void display_matrix() const
    {
        for (const auto &row : matrix)
        {
            for (double elem : row)
            {
                std::cout << elem << " ";
            }
            std::cout << std::endl;
        }
    }
    void initial_multiplicity(std::vector<double> vec)
    {
        multiplicity = vec;
        for (auto &i : multiplicity)
        {
            if (i < 0)
            {
                i = 0;
            }
        }
        normalize(multiplicity);
    }
    std::vector<std::vector<double>> matrix;
    std::vector<double> multiplicity;
    int m, n;
    int N;
    double epsilon0;
};
/// //////////////////////////////////////////////////
/// //////////////////////////////////////////////////
vector<double> multiply_matrix_vector(vector<vector<double>> A, vector<double> v)
{
    vector<double> w(A.size(), 0); // Initialize result vector

    for (size_t i = 0; i < A.size(); ++i)
    {
        for (size_t j = 0; j < v.size(); ++j)
        {
            w[i] += A[i][j] * v[j];
        }
    }

    return w;
}

std::vector<double> calculate_fold(Setup S)
{
    vector<double> out = multiply_matrix_vector(S.matrix, S.multiplicity);
    normalize(out);
    return out;
}

double Loss(std::vector<double> calculated_fold, std::vector<double> measured_fold, std::vector<double> uncertainty)
{
    double out = 0;
    for (int i = 0; i < calculated_fold.size(); i++)
    {
        out += pow((calculated_fold[i] - measured_fold[i]) / uncertainty[i], 2);
    }
    return out;
}
/////////////////////////////
/////////brute force//////////
/////////////////////////////
/// //////////////////////////////
vector<double> grad_success_brute(Setup S, vector<double> measured_fold, vector<double> uncertainty, double dx = 0.001)
{
    vector<double> grad;
    vector<double> copy = S.multiplicity;
    for (size_t i = 0; i < copy.size(); ++i)
    {
        S.initial_multiplicity(copy);
        double x = Loss(calculate_fold(S), measured_fold, uncertainty);
        copy[i] += dx;
        S.initial_multiplicity(copy);
        double new_x = Loss(calculate_fold(S), measured_fold, uncertainty);
        grad.push_back((new_x - x) / dx); // Correct gradient calculation
        copy[i] -= dx;                    // Reset the change
    }
    // print_distribution(grad);
    return grad;
}
void iterate_brute(Setup &S, vector<double> measured_fold, vector<double> uncertainty, double d0 = 0.001, double lr = 0.00001)
{
    vector<double> grad = grad_success_brute(S, measured_fold, uncertainty, d0);
    product_vector_scalar(grad, -lr); // Gradient descent step
    vector<double> copy = S.multiplicity;
    sum_vectors(grad, copy);
    S.initial_multiplicity(copy);
}
std::pair<double, std::vector<double>> optimize(Setup &setup, std::vector<double> measured_fold, std::vector<double> uncertainty, double lr = 10000)
{
    int i = 0;
    double buff = Loss(calculate_fold(setup), measured_fold, uncertainty);
    vector<double> multiplicity0 = setup.multiplicity;
    while (1)
    {
        iterate_brute(setup, measured_fold, uncertainty, 0.001, lr);
        double el = Loss(calculate_fold(setup), measured_fold, uncertainty);
        if (el > buff || isnan(el) or i == 100)
        {
            lr /= 2;
            if (lr < 1e-30)
            {
                break;
            }
            setup.initial_multiplicity(multiplicity0);
            i = 0; // Reset iteration counter to restart optimization with new learning rate
        }
        else
        {
            i++;
            multiplicity0 = setup.multiplicity;
            buff = el;
            // cout << "Loss: " << buff << endl;
        }
    }
    return std::make_pair(buff, multiplicity0);
}
////////////////////////////////////////////////////////////////
std::pair<double, std::vector<double>> parallel_optimize(Setup setup, const std::vector<double> &measured_fold, const std::vector<double> &uncertainty, const std::vector<double> &el, double lr = 10000)
{
    setup.initial_multiplicity(el);
    return optimize(setup, measured_fold, uncertainty, lr);
}

std::pair<double, std::vector<double>> optimum_auxiliary(Setup setup, int grid_shots, const std::vector<double> &measured_fold, const std::vector<double> &uncertainty, double lr = 10000)
{
    std::vector<std::vector<double>> combinations = getCombinations(setup.m, grid_shots);
    double buf = std::numeric_limits<double>::max();
    std::vector<double> out;

    // Use ROOT's TThreadExecutor for parallel execution
    ROOT::TThreadExecutor executor;

    // Define the lambda function to be executed in parallel
    auto task = [&](const std::vector<double> &el)
    {
        return parallel_optimize(setup, measured_fold, uncertainty, el);
    };

    // Execute the task in parallel over all combinations
    auto results = executor.Map(task, combinations);

    // Collect results
    for (const auto &result : results)
    {
        if (result.first < buf)
        {
            buf = result.first;
            out = result.second;
        }
    }
    return std::make_pair(buf, out);
}
pair<double, vector<double>> find_optimum(Setup setup, int grid_shots, vector<double> measured_fold, vector<double> uncertainty, double lr = 10000, double de = 0, double ef_shots = 1)
{
    if (de == 0)
    {
        return optimum_auxiliary(setup, grid_shots, measured_fold, uncertainty);
    }
    double epsilon = setup.epsilon0;
    double buf = __DBL_MAX__;
    vector<double> out;
    for (double ep2 = epsilon - de; ep2 < epsilon + de; ep2 += 2 * de / ef_shots)
    {
        setup.epsilon0 = ep2;
        auto result = optimum_auxiliary(setup, grid_shots, measured_fold, uncertainty);
        if (result.first < buf)
        {
            buf = result.first;
            out = result.second;
        }
    }
    return make_pair(buf, out);
}
pair<vector<double>, vector<double>> fold_distribution(vector<double> fold, int k)
{
    vector<double> result;
    int foldSize = fold.size();
    if (k < foldSize)
    {
        result.assign(fold.begin(), fold.begin() + k);
    }
    if (k >= foldSize)
    {
        result = fold;
        result.resize(k, 0.0);
    }
    vector<double> uncertainties;
    double sum = 0;
    for (auto el : fold)
    {
        sum += el;
    }

    for (int i = 0; i < k; i++)
    {
        if (result[i] == 0)
        {
            uncertainties.push_back(1 / sum);
        }
        else
        {
            uncertainties.push_back(sqrt(result[i]) / sum);
        }
    }

    normalize(result);
    normalize(uncertainties);
    return make_pair(result, uncertainties);
}