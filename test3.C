#include "functions.h"
#include <iostream>
using namespace std;

void test3()
{
    Setup setup;
    int m = 4;             // Maximum multiplicity
    int n = 4;             // Maximum fold value
    int N = 21;            // Example value for N
    double epsilon = 0.11; // Example value for total efficiency
    setup.set_limits(m, n);
    setup.set_epsilon0(epsilon);
    setup.set_N(N);
    setup.fill_matrix();
    setup.initial_multiplicity({1, 7, 7, 3});
    vector<double> real_multiplicity = setup.multiplicity;
    vector<double> measured_fold = calculate_fold(setup);
    vector<double> uncertainty = {0.001, 0.01, 0.1, 1};
    setup.set_epsilon0(epsilon);
    auto result = find_optimum(setup, 9, measured_fold, uncertainty, 0.0, 10);
    setup.initial_multiplicity(result.second);
    cout << "Final response" << endl;
    cout << "Measured fold distribution:" << endl;
    print_distribution(calculate_fold(setup));
    cout << "Calculated fold distribution" << endl;
    print_distribution(measured_fold);
    cout << "Real Multiplicity" << endl;
    print_distribution(real_multiplicity);
    cout << "Calculated multiplicity" << endl;
    print_distribution(setup.multiplicity);
}
