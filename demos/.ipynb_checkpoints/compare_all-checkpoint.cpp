#include "timestepper.hpp"
#include "mass_spring_function.hpp"
#include "ExplicitEuler.hpp"   // <--- ADDED THIS (was missing in your original)
#include "ImprovedEuler.hpp"
#include "ImplicitEuler.hpp"
#include "CrankNicolson.hpp"
#include <memory>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>              // <--- ADDED for std::to_string

// IMPORTANT: If your classes (Vector, MassSpring, etc.) are in a namespace, 
// you must uncomment the line below. Check your header files to see if they use one.
// using namespace asc_ode; 

int main() {
    // 1. Setup the problem (Mass Spring System)
    auto rhs = std::make_shared<MassSpring>(1.0, 1.0);
    
    // 2. Initial Condition
    Vector<> y0{1.0, 0.0};

    // 3. Define the helper lambda to run simulations
    auto run = [&](auto& stepper, const std::string& name, double tau) {
        Vector<> y = y0;
        
        // Construct filename (ensure "results" folder exists!)
        std::string filename = "results/" + name + "_tau" + std::to_string(tau) + ".txt";
        std::ofstream f(filename);
        
        if (!f.is_open()) {
            std::cerr << "Error: Could not open file " << filename << ". Did you create the 'results' folder?\n";
            return;
        }

        f << std::fixed << std::setprecision(10);
        f << "0 " << y(0) << " " << y(1) << "\n";
        
        int steps = int(40.0 / tau + 0.5);
        for (int i = 0; i < steps; ++i) {
            stepper->doStep(tau, y);
            f << (i+1)*tau << " " << y(0) << " " << y(1) << "\n";
        }
        std::cout << name << " t=" << tau << " done\n";
    };

    // 4. Register all solvers
    // Note: ExplicitEuler was causing an error because it wasn't included at the top.
    std::vector<std::pair<std::unique_ptr<TimeStepper>, std::string>> methods;
    methods.push_back({std::make_unique<ExplicitEuler>(rhs), "ExplicitEuler"});
    methods.push_back({std::make_unique<ImprovedEuler>(rhs), "ImprovedEuler"});
    methods.push_back({std::make_unique<ImplicitEuler>(rhs), "ImplicitEuler"});
    methods.push_back({std::make_unique<CrankNicolson>(rhs), "CrankNicolson"});

    // 5. Run simulations for different time steps
    for (double tau : {0.1, 0.5, 1.5, 3.0}) {
        // We have to iterate carefully because unique_ptr cannot be copied, only moved.
        // However, since we stored them in the vector, we can access them by reference.
        for (auto& p : methods) {
            // Reset the stepper if necessary (depends on implementation, usually safe for simple steppers)
            run(p.first, p.second, tau);
        }
    }

    std::cout << "\nALL DONE! Results in results/ folder\n";
}