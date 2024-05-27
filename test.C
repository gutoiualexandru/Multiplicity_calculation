#include "functions.h"
// #include "generators.h"
#include <iostream>
using namespace std;

void test()
{
    Setup setup;
    int m = 4;             // Maximum multiplicity
    int n = 4;             // Maximum fold value
    int N = 21;            // Example value for N
    double epsilon = 0.11; // Example value for total efficiency
    setup.set_limits(m, n);
    setup.set_epsilon0(epsilon + 0.001);
    setup.set_N(N);
    setup.fill_matrix();
    setup.initial_multiplicity({0.1, 7, 7, 3});
    vector<double> real_multiplicity = setup.multiplicity;
    vector<double> measured_fold = calculate_fold(setup);
    vector<double> uncertainty = {0.001, 0.01, 0.1, 1};
    setup.set_epsilon0(epsilon);
    double de = 0.01;
    double samples_e = 100;
    std::vector<std::vector<double>> combinations = getCombinations(m, 6);
    std::unordered_map<double, std::vector<double>> hashmap2;
    for (double ep2 = epsilon - de; ep2 < epsilon + de; ep2 += 2 * de / samples_e)
    {
        setup.set_epsilon0(ep2);
        setup.fill_matrix();
        std::unordered_map<double, std::vector<double>> hashmap;
        for (auto el : combinations)
        {
            setup.initial_multiplicity(el);
            pair<double, std::vector<double>> result = optimize(setup, measured_fold, uncertainty);
            hashmap[result.first] = result.second;
        }
        double ll = findLowestKeyVector(hashmap);
        setup.initial_multiplicity(hashmap[ll]);
        hashmap2[ll] = setup.multiplicity;
    }
    setup.initial_multiplicity(hashmap2[findLowestKeyVector(hashmap2)]);
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
