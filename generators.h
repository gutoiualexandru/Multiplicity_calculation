#include <iostream>
#include <unordered_map>
#include <vector>
#include <limits>

// Helper function to find all combinations
void findCombinations(int m, int sum, int currentSum, std::vector<int> &currentCombination, std::vector<std::vector<int>> &results)
{
    if (currentCombination.size() == m)
    {
        if (currentSum == sum)
        {
            results.push_back(currentCombination);
        }
        return;
    }

    for (int i = 0; i <= sum; ++i)
    {
        if (currentSum + i <= sum)
        {
            currentCombination.push_back(i);
            findCombinations(m, sum, currentSum + i, currentCombination, results);
            currentCombination.pop_back();
        }
    }
}

// Function to get all combinations of m integers summing to sum
std::vector<std::vector<double>> getCombinations(int m, int sum)
{
    std::vector<std::vector<int>> intResults;
    std::vector<int> currentCombination;

    findCombinations(m, sum, 0, currentCombination, intResults);

    std::vector<std::vector<double>> doubleResults;
    for (const auto &combination : intResults)
    {
        std::vector<double> doubleCombination(combination.begin(), combination.end());
        doubleResults.push_back(doubleCombination);
    }

    return doubleResults;
}
// Function to find the vector associated with the lowest key
double findLowestKeyVector(const std::unordered_map<double, std::vector<double>> &hashmap)
{
    double lowestKey = std::numeric_limits<double>::max();
    std::vector<double> result;

    for (const auto &pair : hashmap)
    {
        if (pair.first < lowestKey)
        {
            lowestKey = pair.first;
            result = pair.second;
        }
    }

    return lowestKey;
}
/////////////////////////////////////
void plotHistogram(std::unordered_map<double, std::vector<double>> hashmap)
{
    // Initialize a ROOT application

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Histogram Canvas", 800, 600);

    // Determine the maximum size of vectors
    int max_size = 0;
    for (const auto &pair : hashmap)
    {
        if (pair.second.size() > max_size)
        {
            max_size = pair.second.size();
        }
    }

    // Create a TH2D histogram
    TH2D *hist = new TH2D("hist", "TH2D Histogram;Element Number;Value", max_size, 0, max_size, 100, 0, 1); // Adjust y-axis range as needed

    // Fill the histogram
    for (const auto &pair : hashmap)
    {
        const std::vector<double> &vec = pair.second;
        for (size_t i = 0; i < vec.size(); ++i)
        {
            hist->Fill(i, vec[i]);
        }
    }

    // Draw the histogram
    hist->Draw("COLZ");

    // Update and display the canvas
    c1->Update();
    c1->Draw();
}