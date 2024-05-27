#include "functions.h"
#include <iostream>
using namespace std;

void test4()
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

    for (auto &el : measured_fold)
    {
        el *= 10000;
    }
    print_distribution(measured_fold);
    auto pairs = fold_distribution(measured_fold, measured_fold.size());
    vector<double> uncertainty = pairs.second;
    measured_fold = pairs.first;
    setup.set_epsilon0(epsilon);
    auto result = find_optimum(setup, 1, measured_fold, uncertainty, 0.0, 10);
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
