#include "functions.h"

void test_particle_efficiency()
{
    // Initialize the target distribution
    vector<double> measured_fold = {3.3857e+09, 2.46312e+09, 5.91864e+08, 1.54377e+08, 3.36631e+07, 6.45244e+06, 1.16963e+06, 212373, 40200, 9060, 8769};

    vector<double> target_y = normalize_distribution(measured_fold);

    const size_t dim = target_y.size();
    const size_t num_particles = 1000;
    const size_t num_iterations = 1000;
    const double omega = 0.005;
    const double c1 = 0.5;
    const double c2 = 0.5;
    TFile *tf = new TFile("Optimization_ISOMIX2.root", "RECREATE");
    TTree *tree = new TTree("tree", "tree");

    Setup setup;
    int m = measured_fold.size(); // Maximum multiplicity
    int n = measured_fold.size(); // Maximum fold value
    int N = 21;                   // Example value for N
    setup.set_limits(m, n);
    setup.set_N(N);
    cout << "Real Fold:" << endl;
    print_distribution(target_y);
    double KL, ef;
    vector<double> estimated_fold;
    vector<double> calculated_multiplicity;
    
    tree->Branch("KL", &KL);
    tree->Branch("Efficiency", &ef);
    tree->Branch("Calculated_Fold", &estimated_fold);
    tree->Branch("Calculated_Multiplicity", &calculated_multiplicity);
    double epsilon = 0.145;       // Example value for total efficiency
    for(double epsilon = 0.12; epsilon<0.17;epsilon+=0.0001)
    {
        
        setup.set_epsilon0(epsilon);
        cout << "efficiency: " << epsilon << endl; 
        ef=epsilon;
        setup.fill_matrix();
        // Initialize the swarm
        vector<Particle> particles;
        for (size_t i = 0; i < num_particles; ++i)
        {
            particles.emplace_back(dim);
        }
        vector<double> global_best_position = particles[0].best_position;
        double global_best_error = kl_divergence(global_best_position, target_y, setup);

        // PSO loop
        for (size_t iter = 0; iter < num_iterations; ++iter)
        {
            for (auto &particle : particles)
            {
                random_device rd;
                mt19937 gen(rd());
                uniform_real_distribution<> dis(0, 1);

                // Update velocity
                for (size_t i = 0; i < dim; ++i)
                {
                    double r1 = dis(gen);
                    double r2 = dis(gen);
                    particle.velocity[i] = omega * particle.velocity[i] +
                                        c1 * r1 * (particle.best_position[i] - particle.position[i]) +
                                        c2 * r2 * (global_best_position[i] - particle.position[i]);
                }
                // Update position
                for (size_t i = 0; i < dim; ++i)
                {
                    particle.position[i] += particle.velocity[i];
                }
                particle.position = normalize_distribution(particle.position);
                // Evaluate fitness
                double error = kl_divergence(particle.position, target_y, setup);
                // Update personal best
                if (error < particle.best_error)
                {
                    particle.best_position = particle.position;
                    particle.best_error = error;
                }

                // Update global best
                if (error < global_best_error)
                {
                    global_best_position = particle.position;
                    global_best_error = error;
                }
            }
        } 
        // Output the best solution
        cout << "Calculated Fold:" << endl;
        setup.initial_multiplicity(global_best_position);
        print_distribution(calculate_fold(setup));
        cout << "Calculated Multiplicity: " << endl;
        for (const auto &val : global_best_position)
        {
            cout << val << " ";
        }
        cout << endl;
        cout << "Best KL divergence: " << global_best_error << endl;
        KL=global_best_error;
        estimated_fold=calculate_fold(setup);
        calculated_multiplicity=setup.multiplicity;
        tree->Fill();
    }
    tf->Write();
}