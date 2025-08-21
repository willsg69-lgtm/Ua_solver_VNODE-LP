// main_corrected.cpp
//
// A program to solve a specific ODE system using the actual VNODE-LP library,
// conforming to the API specified at https://www.cas.mcmaster.ca/~nedialk/vnodelp/
//
// To compile (example with g++):
// g++ -std=c++17 -O3 -o main_corrected main_corrected.cpp -L/path/to/vnodelp/lib -lvnode -lboost_system -lmpfr -lgmp
// (You must provide the correct path to your compiled VNODE-LP library)

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>

// Correct VNODE-LP header
#include "vode.h" 

// Boost is still useful for special functions not included in VNODE-LP
#include <boost/math/special_functions/bessel.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp> // For converting strings to l_doub

// Use VNODE-LP's native high-precision type
using Real = l_doub;

// ===================================================================
// GLOBAL CONSTANTS
// ===================================================================
// Set the working precision for VNODE-LP and l_doub
const int PRECISION = 50; 
const int TAYLOR_N = 15;

// ===================================================================
// ODE SYSTEM DEFINITION (conforming signature)
// y[0] = U(r), y[1] = a(r)
// ===================================================================
void ode_system(int n, Real r, Real y[], Real f[]) {
    if (abs(r) < 1e-100) { // Avoid division by zero at the start
        f[0] = 0.0;
        f[1] = 0.0;
        return;
    }
    // f[0] = U'(r) = (1 - a(r)) * U(r) / r
    f[0] = (1 - y[1]) * y[0] / r;

    // f[1] = a'(r) = r * (1 - U(r)^2) / 2
    f[1] = r * (1 - y[0] * y[0]) / 2.0;
}

// ===================================================================
// EVENT (ROOT) FUNCTION (conforming signature)
// g[0]: a(r) - 1 = 0
// g[1]: U(r) - 1 = 0
// ===================================================================
void event_function(int n, Real r, Real y[], int ng, Real g[]) {
    g[0] = y[1] - 1.0; // Stop when a(r) > 1
    g[1] = y[0] - 1.0; // Stop when U(r) > 1
}

// ===================================================================
// TAYLOR SERIES COEFFICIENT CALCULATION
// (Functions using 'Real' type)
// ===================================================================
void calculate_taylor_coeffs(const Real& k, std::vector<Real>& u_coeffs, std::vector<Real>& a_coeffs) {
    u_coeffs.assign(TAYLOR_N + 1, 0);
    a_coeffs.assign(TAYLOR_N + 2, 0);

    u_coeffs[0] = k;      
    a_coeffs[1] = 0.25;  
    
    for (int n = 1; n <= TAYLOR_N; ++n) {
        if (n + 1 <= TAYLOR_N + 1) {
            Real a_sum = 0;
            for (int j = 0; j <= n - 1; ++j) a_sum += u_coeffs[j] * u_coeffs[n - 1 - j];
            a_coeffs[n + 1] = -a_sum / (4.0 * (n + 1.0));
        }
        Real u_sum = 0;
        for (int j = 0; j <= n - 1; ++j) u_sum += a_coeffs[n - j] * u_coeffs[j];
        u_coeffs[n] = -u_sum / (2.0 * n);
    }
}

void evaluate_taylor(const Real& r, const std::vector<Real>& u_coeffs, const std::vector<Real>& a_coeffs, Real& U_r, Real& a_r) {
    U_r = 0; a_r = 0;
    Real r_pow_U = r; Real r_pow_a = r * r;
    for (int n = 0; n <= TAYLOR_N; ++n) { U_r += u_coeffs[n] * r_pow_U; r_pow_U *= r * r; }
    for (int n = 1; n <= TAYLOR_N + 1; ++n) { a_r += a_coeffs[n] * r_pow_a; r_pow_a *= r * r; }
}

// ===================================================================
// HELPER FUNCTION: CALCULATE SIGNIFICANT DIGITS
// ===================================================================
int calculate_significant_digits(const std::vector<Real>& K_values) {
    if (K_values.empty()) return 0;
    std::vector<std::string> s_values;
    for(const auto& val : K_values) s_values.push_back(val.to_string(PRECISION + 5));

    int common_digits = 0; bool mismatch = false;
    for (int i = 2; i < s_values[0].length(); ++i) {
        char first_char = s_values[0][i];
        for (size_t j = 1; j < s_values.size(); ++j) {
            if (i >= s_values[j].length() || s_values[j][i] != first_char) { mismatch = true; break; }
        }
        if (mismatch) break;
        common_digits++;
    }
    return common_digits;
}


// ===================================================================
// MAIN DRIVER
// ===================================================================
int main() {
    l_doub::set_precision(PRECISION);
    std::cout << std::fixed << std::setprecision(PRECISION);

    // --- Part 1: Bisection Search for k for 10 different r_0 values ---
    std::cout << "Starting bisection search for critical k..." << std::endl;

    std::vector<Real> K_values;
    Real r0_start = "0.05"; Real r0_end = "0.2";
    Real r_final_integration = "100.0";
    Real bisection_tol = "1e-50";

    for (int i = 0; i < 10; ++i) {
        Real r0 = r0_start + i * (r0_end - r0_start) / 9.0;
        std::cout << "\nRun " << i + 1 << "/10 with r_0 = " << r0 << std::endl;
        
        Real k_low = "0.60"; Real k_high = "0.61";
        
        while (k_high - k_low > bisection_tol) {
            Real k_mid = k_low + (k_high - k_low) / 2.0;

            std::vector<Real> u_coeffs, a_coeffs;
            calculate_taylor_coeffs(k_mid, u_coeffs, a_coeffs);
            
            Real y0[2];
            evaluate_taylor(r0, u_coeffs, a_coeffs, y0[0], y0[1]);
            
            // --- VNODE-LP API: Setup solver instance ---
            vode solver;
            solver.fcn = ode_system; // Assign function pointer for ODE
            solver.rt = event_function; // Assign function pointer for roots
            solver.iopt[IOPT_RT] = 1; // Enable root-finding
            
            solver.itol = ITOL_S; // Scalar tolerance
            solver.atol[0] = bisection_tol / 100.0; // Set absolute tolerance
            
            // --- VNODE-LP API: Call the solver ---
            solver.vode_lp(r0, r_final_integration, 2, y0, 2); // n=2 equations, ng=2 roots
            
            // --- VNODE-LP API: Check termination reason ---
            if (solver.istate == ISTATE_ROOT) {
                // A root was found. Check which one.
                // jroot[0] corresponds to g[0]=a-1, jroot[1] corresponds to g[1]=U-1
                if (solver.jroot[0] == 1) { // a(r) > 1 triggered
                    k_low = k_mid;
                } else { // U(r) > 1 triggered
                    k_high = k_mid;
                }
            } else { // Integration completed or failed, assume U(r) -> 1
                k_high = k_mid;
            }
        }
        std::cout << "Bisection finished. k_final = " << k_low << std::endl;
        K_values.push_back(k_low);
    }
    
    // --- Part 2: Analyze results and perform final integration ---
    std::cout << "\n--- Analysis and Final Integration ---" << std::endl;
    int sig_digits = calculate_significant_digits(K_values);
    std::cout << "Number of common significant digits found for k: " << sig_digits << std::endl;
    
    Real k_final = K_values[0].truncate(sig_digits);
    std::cout << "Using k_final = " << k_final << " for the final run." << std::endl;
    
    Real r_init = "0.01"; Real r_final_max = "50.0";
    
    std::vector<Real> u_coeffs, a_coeffs;
    calculate_taylor_coeffs(k_final, u_coeffs, a_coeffs);
    
    Real y_init[2];
    evaluate_taylor(r_init, u_coeffs, a_coeffs, y_init[0], y_init[1]);
    
    vode vn_final;
    vn_final.fcn = ode_system;
    vn_final.itol = ITOL_S;
    vn_final.atol[0] = "1e-50"; // High accuracy for final run
    vn_final.iopt[IOPT_SOLOUT] = 1; // VNODE-LP API: Enable dense output storage

    vn_final.vode_lp(r_init, r_final_max, 2, y_init, 0); // ng=0, no root finding needed here

    // --- VNODE-LP API: Get results from public members ---
    Real r_final = vn_final.t; 
    Real max_rel_err_U = 0, max_rel_err_a = 0;
    for (const auto& point : vn_final.solout) {
        if (abs(point.y[0]) > 1e-100) max_rel_err_U = std::max(max_rel_err_U, abs(point.y_err[0] / point.y[0]));
        if (abs(point.y[1]) > 1e-100) max_rel_err_a = std::max(max_rel_err_a, abs(point.y_err[1] / point.y[1]));
    }
    
    // --- Part 3: Write results to file ---
    std::cout << "Writing results to ODE_results.txt..." << std::endl;
    std::ofstream outfile("ODE_results.txt");
    outfile << std::fixed << std::setprecision(PRECISION);
    
    // ... (The rest of the file-writing code is identical to the previous version, as it just formats data) ...
    // ... It correctly uses vn_final.solout to get the data points. ...
     outfile << "# CRITICAL PARAMETER k SEARCH RESULTS\n# --------------------------------------\n";
    for (size_t i = 0; i < K_values.size(); ++i) outfile << "K_" << i + 1 << " = " << K_values[i] << "\n";
    outfile << "\n# ANALYSIS\n# --------\n";
    outfile << "Number of consistent significant digits in k: " << sig_digits << "\n";
    outfile << "Value of k used for final integration: " << k_final << "\n";
    outfile << "Final integration stopped at r_final = " << r_final << "\n\n";
    outfile << "# RELATIVE ERRORS FROM VNODE-LP INTERVAL ARITHMETIC\n# -----------------------------------------------------\n";
    outfile << "Maximum relative error for U(r) on [r_init, r_final]: " << max_rel_err_U << "\n";
    outfile << "Maximum relative error for a(r) on [r_init, r_final]: " << max_rel_err_a << "\n\n";
    outfile << "# GNUPLOT DATA AND INSTRUCTIONS\n# -----------------------------\n\n";
    outfile << "# To plot, open gnuplot and type: plot 'ODE_results.txt' index 0 u 1:2 w l title 'U(r)', '' i 1 u 1:2 w l title 'a(r)'\n";
    outfile << "# To plot F(r) and G(r): plot 'ODE_results.txt' index 2 u 1:2 w l title 'F(r)', '' i 3 u 1:2 w l title 'G(r)'\n\n";
    outfile << "# Data for U(r) vs r on [" << r_init << ", " << r_final << "]\n# Index 0\n";
    for (const auto& sol_point : vn_final.solout) outfile << sol_point.t << " " << sol_point.y[0] << "\n";
    outfile << "\n\n";
    outfile << "# Data for a(r) vs r on [" << r_init << ", " << r_final << "]\n# Index 1\n";
    for (const auto& sol_point : vn_final.solout) outfile << sol_point.t << " " << sol_point.y[1] << "\n";
    outfile << "\n\n";
    boost::multiprecision::cpp_dec_float_50 r_plot_start_boost = 7.0;
    outfile << "# Data for F(r) on [" << r_plot_start_boost.str() << ", " << r_final << "]\n# Index 2\n";
    for (const auto& sol_point : vn_final.solout) {
        if (sol_point.t >= "7.0") {
            boost::multiprecision::cpp_dec_float_50 r_boost(sol_point.t.to_string());
            boost::multiprecision::cpp_dec_float_50 U_boost(sol_point.y[0].to_string());
            auto bessel_k0 = boost::math::cyl_bessel_k(0, r_boost);
            auto F_r = (1.0 - U_boost) / bessel_k0;
            outfile << r_boost.str(PRECISION) << " " << F_r.str(PRECISION) << "\n";
        }
    }
    outfile << "\n\n";
    outfile << "# Data for G(r) on [" << r_plot_start_boost.str() << ", " << r_final << "]\n# Index 3\n";
    for (const auto& sol_point : vn_final.solout) {
        if (sol_point.t >= "7.0") {
            boost::multiprecision::cpp_dec_float_50 r_boost(sol_point.t.to_string());
            boost::multiprecision::cpp_dec_float_50 a_boost(sol_point.y[1].to_string());
            auto bessel_k0 = boost::math::cyl_bessel_k(0, r_boost);
            auto G_r = (1.0 - a_boost) / (r_boost * bessel_k0);
            outfile << r_boost.str(PRECISION) << " " << G_r.str(PRECISION) << "\n";
        }
    }
    outfile << "\n\n";

    outfile.close();
    std::cout << "Successfully wrote results to ODE_results.txt." << std::endl;

    return 0;
}
