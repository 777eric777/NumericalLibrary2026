# NumericalLibrary

Modular numerical methods library. This repository contains a set of
independent modules implementing common algorithms used in numerical
computations: linear systems, interpolation, approximation, integration,
ordinary differential equation (ODE) solvers and nonlinear equation solvers.

Repository layout

NumericalLibrary
├── include/             # Public headers for each module
│   ├── LinearSystems.h
│   ├── Interpolation.h
│   ├── Approximation.h
|   ├── Integration.h
│   ├── ODE.h
│   └── NonlinearEquations.h
├── src/                 # Implementation files
│   ├── LinearSystems.cpp
│   ├── Interpolation.cpp
│   ├── Approximation.cpp
│   ├── Integration.cpp
│   ├── ODE.cpp
│   └── NonlinearEquations.cpp
├── tests/               # unit tests and test runner
├── examples/            # example programs demonstrating library usage
└── README.md			 # [You're here!]

Core modules and notable APIs

- LinearSystems
  - Solvers for linear systems (Gaussian elimination, LU, ...).

- Interpolation
  - Polynomial and piecewise interpolation helpers.

- Approximation
  - Least-squares and regression utilities.

- Integration
  - Numerical quadrature rules (trapezoid, Simpson, adaptive routines).

- ODE
  - Ordinary differential equation solvers (e.g., RK4 implemented as
	rungeKutta4).

- NonlinearEquations
  - Root-finding routines for scalar equations. Current public functions:
	- newton_analytical(df, f, x0, tol, max_iter, outFile)
	  Newton's method when analytical derivative is available.
	- newton_numerical(f, x0, tol, max_iter, outFile)
	  Newton's method using a central-difference approximation of the derivative.
	- secantMethod(f, x0, x1, tolerance, maxIterations, outFile)
	  Secant method requiring two initial approximations.
	- derivativeNumerical(f, x, h)
	  Central-difference numerical derivative helper.
	- bisection(f, a, b, root, tol, maxIter, verbose, iter_log)
	  Robust bracketed root-finder using interval halving.
	- false_position(f, a, b, root, tol, maxIter, verbose, iter_log)
	  Regula falsi (false position) method — bracketed and often faster than bisection.
	- find_roots_bisection / find_roots_false_position
	  Helper routines that subdivide an interval and apply bisection or
	  false-position on subintervals where a sign change (or near-zero sample)
	  is detected. They return a list of distinct roots and optional per-root
	  iteration logs used for diagnostics.

Using the library (quick example)

1. Include the aggregated header in your program:

   #include "NumericalLibrary.h"
   using namespace NumericalLibrary;

2. Call functions directly. Example (Newton analytic):

   auto f = [](double x){ return x*x - 2.0; };
   auto df = [](double x){ return 2.0*x; };
   double root = newton_analytical(df, f, 1.0, 1e-12, 50, "");

Build instructions

The project is compatible with CMake + Ninja and MSVC. Example:

  mkdir build
  cd build
  cmake -G "Ninja" -DCMAKE_BUILD_TYPE=Debug ..
  cmake --build . --config Debug

Alternatively open the repository in Visual Studio (File > Open > CMake)
and build the desired target.

Notes regarding testing:

- Tests: The tests/ folder contains two diagnostical tests for each of the functions in NumericalLibrary.
  Feel free to explore the edge cases and the general application of those methods.
- Examples: examples/ contains usage_example.cpp which demonstrates basic
  Newton and RK4 usage, as well as demo_odes.cpp for rangeKutta4 and eulerExplicit and demo_root_finding.cpp
  with Newton (with two types of derivatives) and the secantMethod

