// main.cpp
//
// A program to solve a specific ODE system using VNODE-LP,
// performing a bisection search for a critical parameter k.
//
// To compile (example with g++):
// g++ -std=c++17 -O3 -o main main.cpp -lvnode -lboost_system -lboost_thread -lmpfr -lgmp
// (The exact linking flags for VNODE-LP and Boost may vary based on your system setup)

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>

// VNODE-LP header. The actual header name might differ.
#include <vnode.h> 

// Boost libraries for high-precision floating point numbers and special functions.
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/bessel.hpp>

// Use a 50-digit decimal float type from Boost. 
// VNODE-LP might have its own high-precision type (e.g., l_doub),
// in which case this should be replaced with that type.
using mp_float = boost::multiprecision::cpp_dec_float_50;

// ===================================================================
// GLOBAL CONSTANTS
// ===================================================================
const int PRECISION = 50;
const int TAYLOR_N = 15; // N=15 means up to r^31 for U and r^32 for a.

// ===================================================================
// ODE SYSTEM DEFINITION
// y[0] = U(r), y[1] = a(r)
// ===================================================================
void ode_system(int n, const mp_float& r, const mp_float y[], mp_float f[]) {
    // f[0] = U'(r) = (1 - a(r)) * U(r) / r
    f[0] = (1 - y[1]) * y[0] / r;

    // f[1] = a'(r) = r * (1 - U(r)^2) / 2
    f[1] = r * (1 - y[0] * y[0]) / 2.0;
}

// ===================================================================
// EVENT FUNCTIONS FOR SOLVER TERMINATION
// The solver stops when any of these functions cross zero.
// g[0]: a(r) - 1 = 0
// g[1]: U(r) - 1 = 0
// ===================================================================
void event_function(int n, const mp_float& r, const mp_float y[], int n_events, mp_float g[]) {
    g[0] = y[1] - 1.0; // Stop when a(r) > 1
    g[1] = y[0] - 1.0; // Stop when U(r) > 1
}

// ===================================================================
// TAYLOR SERIES COEFFICIENT CALCULATION
// ===================================================================
void calculate_taylor_coeffs(const mp_float& k,
                             std::vector<mp_float>& u_coeffs,
                             std::vector<mp_float>& a_coeffs) {
    u_coeffs.assign(TAYLOR_N + 1, 0);
    a_coeffs.assign(TAYLOR_N + 2, 0);

    // Initial conditions from the problem statement
    u_coeffs[0] = k;      // Corresponds to u_1
    a_coeffs[1] = 0.25;   // Corresponds to a_2

    // Iteratively compute coefficients
    for (int n = 1; n <= TAYLOR_N; ++n) {
        // Compute a_{2(n+1)}
        if (n + 1 <= TAYLOR_N + 1) {
            mp_float a_sum = 0;
            // Summation runs from j=0 to (n+1)-2 = n-1
            for (int j = 0; j <= n - 1; ++j) {
                a_sum += u_coeffs[j] * u_coeffs[n - 1 - j];
            }
            a_coeffs[n + 1] = -a_sum / (4.0 * (n + 1.0));
        }

        // Compute u_{2n+1}
        mp_float u_sum = 0;
        // Summation runs from j=0 to n-1
        for (int j = 0; j <= n - 1; ++j) {
            u_sum += a_coeffs[n - j] * u_coeffs[j];
        }
        u_coeffs[n] = -u_sum / (2.0 * n);
    }
}


// ===================================================================
// TAYLOR SERIES EVALUATION
// ===================================================================
void evaluate_taylor(const mp_float& r,
                     const std::vector<mp_float>& u_coeffs,
                     const std::vector<mp_float>& a_coeffs,
                     mp_float& U_r, mp_float& a_r) {
    U_r = 0;
    a_r = 0;
    mp_float r_pow_U = r;
    mp_float r_pow_a = r * r;

    // Horner-like scheme is more stable, but direct summation is fine for small r
    for (int n = 0; n <= TAYLOR_N; ++n) {
        U_r += u_coeffs[n] * r_pow_U;
        r_pow_U *= r * r;
    }

    for (int n = 1; n <= TAYLOR_N + 1; ++n) {
        a_r += a_coeffs[n] * r_pow_a;
        r_pow_a *= r * r;
    }
}


// ===================================================================
// HELPER FUNCTION: CALCULATE SIGNIFICANT DIGITS
// ===================================================================
int calculate_significant_digits(const std::vector<mp_float>& K_values) {
    if (K_values.empty()) return 0;
    
    std::vector<std::string> s_values;
    for(const auto& val : K_values) {
        s_values.push_back(val.str(PRECISION + 5, std::ios_base::fixed));
    }

    int common_digits = 0;
    bool mismatch = false;
    // Start after "0."
    for (int i = 2; i < s_values[0].length(); ++i) {
        char first_char = s_values[0][i];
        for (size_t j = 1; j < s_values.size(); ++j) {
            if (i >= s_values[j].length() || s_values[j][i] != first_char) {
                mismatch = true;
                break;
            }
        }
        if (mismatch) {
            break;
        }
        common_digits++;
    }
    return common_digits;
}


// ===================================================================
// MAIN DRIVER
// ===================================================================
int main() {
    std::cout << std::fixed << std::setprecision(PRECISION);

    // --- Part 1: Bisection Search for k for 10 different r_0 values ---
    std::cout << "Starting bisection search for critical k..." << std::endl;

    std::vector<mp_float> K_values;
    mp_float r0_start = "0.05";
    mp_float r0_end = "0.2";
    mp_float r_final_integration = "100.0";
    mp_float bisection_tol = "1e-50";

    for (int i = 0; i < 10; ++i) {
        mp_float r0 = r0_start + i * (r0_end - r0_start) / 9.0;
        std::cout << "\nRun " << i + 1 << "/10 with r_0 = " << r0 << std::endl;
        
        mp_float k_low = "0.60";
        mp_float k_high = "0.61";
        
        int p = 0;
        while (k_high - k_low > bisection_tol) {
            mp_float k_mid = k_low + (k_high - k_low) / 2.0;

            // 1. Get initial conditions at r0 for this k_mid
            std::vector<mp_float> u_coeffs, a_coeffs;
            calculate_taylor_coeffs(k_mid, u_coeffs, a_coeffs);
            
            mp_float U0, a0;
            evaluate_taylor(r0, u_coeffs, a_coeffs, U0, a0);
            
            mp_float y0[2] = {U0, a0};

            // 2. Setup VNODE-LP solver
            vnode<mp_float> vn;
            vn.set_function(ode_system);
            vn.set_event_function(2, event_function); // 2 events
            
            // Set high precision tolerances
            vn.set_tolerances(bisection_tol / 100.0, bisection_tol / 100.0);
            
            // 3. Integrate
            // The actual VNODE-LP call might look different.
            // This is a representative syntax.
            // We expect it to return a status indicating which event triggered.
            int term_reason = vn.integrate(r_final_integration, r0, 2, y0);
            
            // 4. Update k interval based on termination reason
            // We assume vn.get_event_info() or similar tells us which event happened.
            // Let's assume term_reason == 1 for event 0 (a>1) and term_reason == 2 for event 1 (U>1)
            // This logic needs to be adapted to the actual VNODE-LP API return codes.
            if (term_reason == 1) { // a(r) > 1 triggered
                k_low = k_mid;
            } else { // U(r) > 1 triggered or r_final_integration reached
                k_high = k_mid;
            }
            p++;
        }
        std::cout << "Bisection finished in " << p << " steps. k_final = " << k_low << std::endl;
        K_values.push_back(k_low);
    }
    
    // --- Part 2: Analyze results and perform final integration ---
    std::cout << "\n--- Analysis and Final Integration ---" << std::endl;
    
    int sig_digits = calculate_significant_digits(K_values);
    std::cout << "Number of common significant digits found for k: " << sig_digits << std::endl;
    
    // Truncate K_1 to the number of significant digits
    std::string k1_str = K_values[0].str(PRECISION + 5, std::ios_base::fixed);
    std::string k_final_str = k1_str.substr(0, 2 + sig_digits);
    mp_float k_final = k_final_str;
    
    std::cout << "Using k_final = " << k_final << " for the final run." << std::endl;
    
    // --- Final integration run ---
    mp_float r_init = "0.01";
    mp_float r_final_max = "50.0";
    
    std::vector<mp_float> u_coeffs, a_coeffs;
    calculate_taylor_coeffs(k_final, u_coeffs, a_coeffs);
    
    mp_float U_init, a_init;
    evaluate_taylor(r_init, u_coeffs, a_coeffs, U_init, a_init);
    
    mp_float y_init[2] = {U_init, a_init};
    
    vnode<mp_float> vn_final;
    vn_final.set_function(ode_system);
    vn_final.set_event_function(2, event_function);
    
    // Request solution storage and interval arithmetic for error bounds
    vn_final.set_solution_output(SOL_STORE); // Assuming this enum/flag exists
    vn_final.enable_interval_arithmetic();   // Assuming such a method exists

    vn_final.integrate(r_final_max, r_init, 2, y_init);
    
    // Get final results
    mp_float r_final = vn_final.get_final_time(); // Get actual stopping time
    auto solution = vn_final.get_solution(); // std::vector<std::pair<mp_float, mp_float*>>
    auto errors = vn_final.get_error_bounds(); // Assuming returns error bounds
    
    // Estimate relative errors. This is highly dependent on VNODE-LP's API.
    // Here we simulate by finding the max relative error bound over the solution.
    mp_float max_rel_err_U = 0;
    mp_float max_rel_err_a = 0;
    for (size_t i = 0; i < solution.size(); ++i) {
        mp_float U_val = solution[i].second[0];
        mp_float a_val = solution[i].second[1];
        mp_float U_err = errors[i].second[0]; // Width of interval for U
        mp_float a_err = errors[i].second[1]; // Width of interval for a
        
        if (abs(U_val) > 1e-50) {
            max_rel_err_U = std::max(max_rel_err_U, abs(U_err / U_val));
        }
        if (abs(a_val) > 1e-50) {
            max_rel_err_a = std::max(max_rel_err_a, abs(a_err / a_val));
        }
    }


    // --- Part 3: Write results to file ---
    std::cout << "Writing results to ODE_results.txt..." << std::endl;
    std::ofstream outfile("ODE_results.txt");
    outfile << std::fixed << std::setprecision(PRECISION);
    
    outfile << "# CRITICAL PARAMETER k SEARCH RESULTS\n";
    outfile << "# --------------------------------------\n";
    for (size_t i = 0; i < K_values.size(); ++i) {
        outfile << "K_" << i + 1 << " = " << K_values[i] << "\n";
    }
    outfile << "\n";
    
    outfile << "# ANALYSIS\n";
    outfile << "# --------\n";
    outfile << "Number of consistent significant digits in k: " << sig_digits << "\n";
    outfile << "Value of k used for final integration: " << k_final << "\n";
    outfile << "Final integration stopped at r_final = " << r_final << "\n\n";

    outfile << "# RELATIVE ERRORS FROM VNODE-LP INTERVAL ARITHMETIC\n";
    outfile << "# -----------------------------------------------------\n";
    outfile << "Maximum relative error for U(r) on [r_init, r_final]: " << max_rel_err_U << "\n";
    outfile << "Maximum relative error for a(r) on [r_init, r_final]: " << max_rel_err_a << "\n\n";


    outfile << "# GNUPLOT DATA AND INSTRUCTIONS\n";
    outfile << "# -----------------------------\n\n";
    outfile << "# To plot, open gnuplot and type: plot 'ODE_results.txt' index 0 u 1:2 w l title 'U(r)', '' i 1 u 1:2 w l title 'a(r)'\n";
    outfile << "# To plot F(r) and G(r): plot 'ODE_results.txt' index 2 u 1:2 w l title 'F(r)', '' i 3 u 1:2 w l title 'G(r)'\n\n";

    // --- Data for U(r) ---
    outfile << "# Data for U(r) vs r on [" << r_init << ", " << r_final << "]\n";
    outfile << "# Index 0\n";
    for (const auto& sol_point : solution) {
        outfile << sol_point.first << " " << sol_point.second[0] << "\n";
    }
    outfile << "\n\n"; // End of data block for gnuplot

    // --- Data for a(r) ---
    outfile << "# Data for a(r) vs r on [" << r_init << ", " << r_final << "]\n";
    outfile << "# Index 1\n";
    for (const auto& sol_point : solution) {
        outfile << sol_point.first << " " << sol_point.second[1] << "\n";
    }
    outfile << "\n\n";

    // --- Data for F(r) and G(r) ---
    mp_float r_plot_start = 7.0;

    // F(r)
    outfile << "# Data for F(r) = (1-U(r))/BesselK(0,r) on [" << std::max(r_init, r_plot_start) << ", " << r_final << "]\n";
    outfile << "# Index 2\n";
     for (const auto& sol_point : solution) {
        if (sol_point.first >= r_plot_start) {
            mp_float r = sol_point.first;
            mp_float U = sol_point.second[0];
            mp_float bessel_k0 = boost::math::cyl_bessel_k(0, r);
            mp_float F_r = (1.0 - U) / bessel_k0;
            outfile << r << " " << F_r << "\n";
        }
    }
    outfile << "\n\n";

    // G(r)
    outfile << "# Data for G(r) = (1-a(r))/(r*BesselK(0,r)) on [" << std::max(r_init, r_plot_start) << ", " << r_final << "]\n";
    outfile << "# Index 3\n";
    for (const auto& sol_point : solution) {
        if (sol_point.first >= r_plot_start) {
            mp_float r = sol_point.first;
            mp_float a = sol_point.second[1];
            mp_float bessel_k0 = boost::math::cyl_bessel_k(0, r);
            mp_float G_r = (1.0 - a) / (r * bessel_k0);
            outfile << r << " " << G_r << "\n";
        }
    }
    outfile << "\n\n";

    outfile.close();
    std::cout << "Successfully wrote results to ODE_results.txt." << std::endl;

    return 0;
}
