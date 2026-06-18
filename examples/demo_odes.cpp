#include "NumericalLibrary.h"
#include <iostream>

// Demo: ODE solvers
// This example demonstrates use of rungeKutta4 and eulerExplicit to integrate
// the exponential decay ODE: y' = -k * y with y(0)=1. Analytical solution: y(t)=exp(-k t)
// Arguments:
// - y0: initial condition
// - t0: initial time
// - t_end: final time
// - h: time step size
// - f: right-hand-side function f(t,y)

int main() {
	using namespace NumericalLibrary;

	double k = 1.0; // decay rate
	auto rhs = [k](double t, double y) { return -k * y; };

	double y0 = 1.0;
	double t0 = 0.0;
	double t_end = 2.0;
	double h = 0.1;

	auto [t_rk, y_rk] = rungeKutta4(y0, t0, t_end, h, rhs);
	std::cout << "RK4 last value y(" << t_rk.back() << ") = " << y_rk.back() << " (analytical=" << exp(-k * t_rk.back()) << ")\n";

	auto [t_eu, y_eu] = eulerExplicit(rhs, t0, y0, h, t_end);
	std::cout << "Euler last value y(" << t_eu.back() << ") = " << y_eu.back() << "\n";

	return 0;
}
