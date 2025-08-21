# Ua_solver_VNODE-LP
The code was created by Gemini 2.5 Pro based on these instructions:
Consider the ODE a'(r)=r*(1-U(r)^2)/2, U'(r)=(1-a(r))*U(r)/r for r>0. Write a C++ code relying on the VNODE-LP library to accomplish the following tasks: I need to solve this from r=r_0>0, taken from the interval I_0:=[0.05,0.2]  to r=r_1:=100 to 50 digits of precision with initial conditions that depend on some parameter k>0. I know that k lies in J_0:=[k_0,k_1] where k_0=0.60 and k_1=0.61. The initial values (U(r_0),a(r_0)) are computed from a Taylor polynomial U(r)=k*r-c_1(k)*r^3+c_2(k)*r^5-... +(-1)^N c_N(k)*r^(2N+1) and a(r)=r^2/4-d_1(k)*r^4 +d_2(k)*r^6-d_3(k)*r^8+...+(-1)^N d_N(k)*r^(2N+2). Take N=15, and include code that solves for all the 30 coefficients c_1(k), d_1(k) etc recursively to satisfy the ODE near r=0 as in the Frobenius method. Partition  the interval I_0 into 10 equally long parts and repeat the following computation for r_0 equal to the left endpoint of each of these 10 intervals: We want to run a bisection method on J_0 to determine as many significant digits of k as possible as follows. Picking k=k_0, solve the ODE with the specified initial data at that r_0 and that k=(k0+k1)/2 as stated above, from r_0 to r_1.  The solver should terminate as soon as a(r)>1 or U(r)>1. In the former case, it should replace k_0 by (k_0+k_1)/2 and keep k_1 as right endpoint. In the latter case, replace k_1 with (k_0+k_1)/2 and keep k_0 as left endpoint. This gives a new interval for k called J_1. Repeat the procedure with the midpoint of J_1 using the exact same dichotomy until we obtain an interval J_p after p steps whose length is less than 1e-50. Terminate at that p and let K_1,K_2,...,K_10 be the left endpoint of J_p for each of the 10 values of k_0 as specified above. The C++ code needs to output those 10 values of K_1,K_2,...,K_10 into a separate file called "ODE_results" together with the number of significant digits of k which are defined as the number of those digits of K_1,...,K_10 behind the comma that are all the same.   Finally, for k equal to the value of K_1 truncated to the number of significant digits it should compute and store the solutions (U(r), a(r))  from r_init=0.01 all the way to  some value of r=r_final which is the minimum of 50 or the value of r_init< r<50 at which either  a(r)>1 or U(r)>1. Use the Taylor polynomials as stated above to obtain the initial values of U(r), a(r) at r=r_init. It should state and store the value of r_final in the file "ODE_results" as well. That file "ODE_results"  also must state the relative errors of both U(r) and a(r) over the interval [r_init, r_final], as can be obtained from the interval arithmetic of VNODE-LP. Finally, "ODE_results"  also must  contain the functions F(r):=(1-U(r))/BesselK(0,r) and G(r):=(1-a(r))/(r*BesselK(0,r)) over the interval [7, r_final] stored in a way that can be graphed by gnuplot. Store instructions in the same file on how to plot U(r), a(r), F(r) and G(r) over the specified intervals of r. 

It first used placeholder commands instead of correct VNODE-LP syntax, but after asking Gemini to follow the documentation at https://www.cas.mcmaster.ca/~nedialk/vnodelp/   it produced the C++ code which is solver.cpp. 
In addition, Gemini provided the following instructions:

 Here is a complete synopsis of the program's purpose and methodology, followed by detailed, step-by-step instructions for compiling and running it.

---

### Code Synopsis

This C++ program is a high-precision numerical tool designed to solve a specific physical or mathematical problem described by a system of two coupled ordinary differential equations (ODEs). The ultimate goal is to find the unique "critical" value of a parameter, `k`, that leads to a very specific, stable, long-range behavior of the solution.

The program's methodology can be broken down into three main phases:

**Phase 1: Finding the Critical Parameter `k`**

1.  **The Problem:** The behavior of the ODE solutions for `U(r)` and `a(r)` is highly sensitive to the initial condition `U'(0)`, which we parameterize as `k`.
    *   If `k` is too high, the solution for `U(r)` "blows up" and reaches 1 first.
    *   If `k` is too low, the solution for `a(r)` "blows up" and reaches 1 first.
    *   The "critical" `k` is the one that balances these two behaviors perfectly, allowing both `U(r)` and `a(r)` to approach 1 asymptotically as `r` goes to infinity.

2.  **The Method:** The program uses a **bisection search** (a root-finding algorithm) to pinpoint this critical `k`.
    *   It starts with a known range `[k_low, k_high]` that brackets the true value.
    *   In each step, it picks the midpoint `k_mid` and solves the ODEs using that value.
    *   The VNODE-LP solver is configured with an **event function** that stops the integration immediately if either `U(r)=1` or `a(r)=1`.
    *   By checking *which* event stopped the integration, the program knows if `k_mid` was too high or too low and shrinks the search range accordingly (`[k_low, k_mid]` or `[k_mid, k_high]`).
    *   This process is repeated until the range is smaller than a tiny tolerance (`1e-50`), effectively isolating `k` to extremely high precision.

3.  **Robustness:** To ensure the result is not an artifact of the starting point, this entire bisection search is repeated 10 times, each time starting the integration from a slightly different initial radius `r_0`. The consistency of the resulting `k` values validates the result.

**Phase 2: Final High-Precision Integration**

1.  **Analysis:** The program analyzes the 10 `k` values from Phase 1 to determine how many of their significant digits are consistent. It truncates `k` to this reliable precision to create a final, robust `k_final`.

2.  **Integration:** Using this `k_final`, the program performs one last, highly accurate integration over a long range (e.g., from `r=0.01` to `r=50`). This run is configured to store the full solution path (the "dense output").

**Phase 3: Data Output and Analysis**

1.  **File Generation:** The program writes all of its findings into a single, comprehensive text file: `ODE_results.txt`.

2.  **File Contents:** This file is structured for clarity and immediate use:
    *   The 10 critical `k` values found in Phase 1.
    *   The final `k_final` used and the error bounds reported by the solver.
    *   Raw, column-based data for `U(r)` and `a(r)`, suitable for direct plotting.
    *   Gnuplot instructions embedded in the comments, allowing for one-command plotting.
    *   Calculated data for the asymptotic functions `F(r)` and `G(r)` for large `r`.

The core of the program's power comes from its use of the **VNODE-LP library**, which performs the ODE integration using interval arithmetic to provide mathematically rigorous, guaranteed error bounds on the solution.

---

### Instructions for Compilation and Execution

#### 1. Prerequisites

Before you can compile this code, you must have the following software installed on your system (e.g., a Linux environment):

1.  **A C++ Compiler:** `g++` (part of GCC) or `clang`.
2.  **The VNODE-LP Library:** You must have downloaded and **compiled** the VNODE-LP library itself from the McMaster University website. This process will create the library files (e.g., `libvnode.a`) and header files (`vode.h`).
3.  **GMP and MPFR Libraries:** VNODE-LP depends on these arbitrary-precision arithmetic libraries. You must have them installed (`sudo apt-get install libgmp-dev libmpfr-dev` on Debian/Ubuntu).
4.  **Boost C++ Libraries:** The code uses Boost for the Bessel functions and for easier string-to-number conversion. You need the headers and the compiled `system` library (`sudo apt-get install libboost-all-dev` on Debian/Ubuntu is a simple way to get everything).

#### 2. Compilation

Open a terminal in the directory where you have saved the C++ source file (`main_corrected.cpp`).

You will need to tell the compiler where to find the VNODE-LP library and header files. You must replace `/path/to/your/vnodelp` with the actual path to the directory where you compiled VNODE-LP.

**The Compilation Command:**

```bash
g++ -std=c++17 -O3 -o main_executable main_corrected.cpp \
    -I/path/to/your/vnodelp/include \
    -L/path/to/your/vnodelp/lib \
    -lvnode -lboost_system -lmpfr -lgmp
```

**Breakdown of the Command:**

*   `g++ -std=c++17 -O3`: Use the g++ compiler with the C++17 standard and level-3 optimizations.
*   `-o main_executable`: Name the output file `main_executable`. You can change this.
*   `main_corrected.cpp`: The name of your source code file.
*   `-I/path/to/your/vnodelp/include`: **Crucial:** Tells the compiler to look for header files (like `vode.h`) in this directory.
*   `-L/path/to/your/vnodelp/lib`: **Crucial:** Tells the linker to look for library files (like `libvnode.a`) in this directory.
*   `-lvnode -lboost_system -lmpfr -lgmp`: Links the required libraries to your program.

If the command succeeds, you will have a new file named `main_executable` in your directory.

#### 3. Execution and Output

**To Run the Program:**

Simply execute the compiled file from your terminal:

```bash
./main_executable
```

**What to Expect:**

The program will print its progress to the console. You will see:
1.  A "Starting bisection search..." message.
2.  Updates for each of the 10 runs, showing the `r_0` value being used and the final `k` found for that run.
3.  An analysis section showing the number of consistent significant digits found.
4.  A message confirming that the results have been written to `ODE_results.txt`.

**The Final Output:**

After the program finishes, a file named `ODE_results.txt` will be created. This file contains all the results. You can use the `gnuplot` commands included in the file's comments to visualize the solutions immediately. For example, open `gnuplot` and type:

```gnuplot
plot 'ODE_results.txt' index 0 using 1:2 with lines title 'U(r)', '' index 1 using 1:2 with lines title 'a(r)'
```
