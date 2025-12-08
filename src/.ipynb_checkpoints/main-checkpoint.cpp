#include <iostream>
#include <fstream>  // For file output
#include <vector>
#include <cmath>    // For math functions like sqrt, fabs, etc.

// Your existing includes
#include "nonlinfunc.hpp"
#include "mass_spring_function.hpp"
#include "newton.hpp"

int main() {
    // Define the system and initial conditions (this might already be in your code)
    std::vector<double> x = {1.0, 0.0};  // Example initial conditions
    double dt = 0.01;  // Time step
    double t_max = 10.0;  // Max time
    int steps = static_cast<int>(t_max / dt);

    // Create a file to output data
    std::ofstream output_file("simulation_data.csv");
    output_file << "Time, x1, x2\n";  // Header for CSV

    // Run the simulation (just an example loop)
    for (int step = 0; step < steps; ++step) {
        // Call the appropriate time-stepping method here (e.g., ExpEuler, ImpEuler)
        // Example of using a simple explicit Euler method:
        std::vector<double> f_val = {1.0, 0.0};  // This should be computed by your model

        for (int i = 0; i < x.size(); ++i) {
            x[i] += dt * f_val[i];  // Update state using Euler method
        }

        // Save the current time and state to the CSV file
        double time = step * dt;
        output_file << time << ", " << x[0] << ", " << x[1] << "\n";
    }

    output_file.close();  // Close the file after writing

    return 0;
}
