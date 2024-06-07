#include "functions.h"
#include <iostream>
using namespace std;

void isomix()
{

    vector<double> measured_fold = {3.3857e+09, 2.46312e+09, 5.91864e+08, 1.54377e+08, 3.36631e+07, 6.45244e+06, 1.16963e+06, 212373, 40200, 9060, 8769};
    // vector<double> measured_fold = {1000.0, 1000.0, 50.0};

    auto pairs = fold_distribution(measured_fold, measured_fold.size());
    vector<double> uncertainty = pairs.second;
    // print_distribution(uncertainty);
    measured_fold = pairs.first;
    Setup setup;
    int m = measured_fold.size(); // Maximum multiplicity
    int n = measured_fold.size(); // Maximum fold value
    int N = 21;                   // Example value for N
    double epsilon = 0.10;        // Example value for total efficiency
    setup.set_limits(m, n);
    setup.set_epsilon0(epsilon);
    setup.set_N(N);
    setup.fill_matrix();
    vector<double> initial = {0, 0, 0, 0, 0, 0.0606124, 0.115636, 0.161354, 0.196866, 0.223154, 0.242377};
    setup.initial_multiplicity(initial);
    auto result = optimize(setup, measured_fold, uncertainty, 0.001, Loss_information);
    // auto result = find_optimum(setup, 5, measured_fold, uncertainty, 0.0001, 0.0, 10, Loss_information);

    setup.initial_multiplicity(result.second);
    cout << "Final response" << endl;
    cout << "Loss:" << result.first << endl;
    cout << "Measured fold distribution:" << endl;
    print_distribution(measured_fold);
    cout << "Calculated fold distribution" << endl;
    print_distribution(calculate_fold(setup));
    cout << "Calculated multiplicity" << endl;
    print_distribution(setup.multiplicity);
}
