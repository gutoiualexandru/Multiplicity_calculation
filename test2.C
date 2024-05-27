#include "generators.h"

void test2()
{
    int m = 3;
    int sum = 5;

    std::vector<std::vector<double>> combinations = getCombinations(m, sum);

    for (const auto &combination : combinations)
    {
        for (double value : combination)
        {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
}
