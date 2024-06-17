#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <random>

using namespace std;
// PSO parameters

// Function to calculate KL divergence
double kl_divergence(const vector<double> &x, const vector<double> &y)
{
    double divergence = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        if (x[i] > 0)
        {
            divergence += x[i] * log(x[i] / y[i]);
        }
    }
    return divergence;
}

// Function to normalize a distribution to ensure it sums to 1
vector<double> normalize_distribution(const vector<double> &dist)
{
    double sum = accumulate(dist.begin(), dist.end(), 0.0);
    vector<double> normalized(dist.size());
    transform(dist.begin(), dist.end(), normalized.begin(),
              [sum](double d)
              { return d / sum; });
    return normalized;
}

// Function to generate a random distribution
vector<double> random_distribution(size_t dim)
{
    vector<double> dist(dim);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);
    generate(dist.begin(), dist.end(), [&]()
             { return dis(gen); });
    return normalize_distribution(dist);
}

// Particle class
class Particle
{
public:
    vector<double> position;
    vector<double> velocity;
    vector<double> best_position;
    double best_error;

    Particle(size_t dim)
        : position(random_distribution(dim)),
          velocity(dim, 0.0),
          best_position(position),
          best_error(numeric_limits<double>::max()) {}
};

// int main() {
//     // Initialize the target distribution
//     vector<double> target_y = random_distribution(dim);

//     // Initialize the swarm
//     vector<Particle> particles;
//     for (size_t i = 0; i < num_particles; ++i) {
//         particles.emplace_back(dim);
//     }

//     vector<double> global_best_position = particles[0].best_position;
//     double global_best_error = kl_divergence(global_best_position, target_y);

//     // PSO loop
//     for (size_t iter = 0; iter < num_iterations; ++iter) {
//         for (auto& particle : particles) {
//             random_device rd;
//             mt19937 gen(rd());
//             uniform_real_distribution<> dis(0, 1);

//             // Update velocity
//             for (size_t i = 0; i < dim; ++i) {
//                 double r1 = dis(gen);
//                 double r2 = dis(gen);
//                 particle.velocity[i] = omega * particle.velocity[i] +
//                                        c1 * r1 * (particle.best_position[i] - particle.position[i]) +
//                                        c2 * r2 * (global_best_position[i] - particle.position[i]);
//             }

//             // Update position
//             for (size_t i = 0; i < dim; ++i) {
//                 particle.position[i] += particle.velocity[i];
//             }
//             particle.position = normalize_distribution(particle.position);

//             // Evaluate fitness
//             double error = kl_divergence(particle.position, target_y);

//             // Update personal best
//             if (error < particle.best_error) {
//                 particle.best_position = particle.position;
//                 particle.best_error = error;
//             }

//             // Update global best
//             if (error < global_best_error) {
//                 global_best_position = particle.position;
//                 global_best_error = error;
//             }
//         }
//     }

//     // Output the best solution
//     cout << "Best position: ";
//     for (const auto& val : global_best_position) {
//         cout << val << " ";
//     }
//     cout << endl;
//     cout << "Best KL divergence: " << global_best_error << endl;

//     return 0;
// }
