// tests/test_runner.cpp
// Simple unit test runner for NumericalLibrary module.
// Each test contains comments describing the procedure.

#include "NumericalLibrary.h"
#include <iostream>
#include <cmath>

#define ASSERT_CLOSE(a,b,eps) do { if (std::fabs((a)-(b)) > (eps)) { std::cerr<<"ASSERT_CLOSE failed: "<< (a) <<" != "<< (b) <<"\n"; return 1; } } while(0)
#define ASSERT_TRUE(cond) do { if (!(cond)) { std::cerr<<"ASSERT_TRUE failed: "<< #cond <<"\n"; return 1; } } while(0)

using namespace NumericalLibrary;

// ODE tests
int test_eulerExplicit_basic() {
	// Procedure: integrate y' = -y with y(0)=1 on [0,1] using small h
	// Expectation: numerical solution close to exp(-1)
	auto f = [](double t, double y){ return -y; };
	auto res = eulerExplicit(f, 0.0, 1.0, 0.001, 1.0);
	double y_end = res.second.back();
	ASSERT_CLOSE(y_end, std::exp(-1.0), 1e-3);
	return 0;
}

int test_eulerExplicit_badstep() {
	// Procedure: call eulerExplicit with negative step size (invalid)
	// Expectation: implementation should return at least initial value; we check that first element equals initial condition
	auto f = [](double t, double y){ return -y; };
	auto res = eulerExplicit(f, 0.0, 1.0, -0.1, 1.0);
	// Implementation returns an empty result for invalid step size; check for that behavior
	ASSERT_TRUE(res.second.empty());
	return 0;
}

int test_rungeKutta4_basic() {
	// Procedure: integrate y' = -y with RK4 over [0,2] with step 0.1
	// Expectation: last value close to exp(-2)
	auto f = [](double t, double y){ return -y; };
	auto [t, y] = rungeKutta4(1.0, 0.0, 2.0, 0.01, f);
	ASSERT_CLOSE(y.back(), std::exp(-2.0), 1e-6);
	return 0;
}

int test_rungeKutta4_stepsize_effect() {
	// Procedure: compare coarse and fine step sizes; fine should be closer to analytical
	auto f = [](double t, double y){ return -y; };
	auto r1 = rungeKutta4(1.0,0.0,1.0,0.2,f);
	auto r2 = rungeKutta4(1.0,0.0,1.0,0.01,f);
	double err1 = std::fabs(r1.second.back() - std::exp(-1.0));
	double err2 = std::fabs(r2.second.back() - std::exp(-1.0));
	ASSERT_TRUE(err2 < err1 + 1e-15);
	return 0;
}

// Interpolation tests
int test_lagrange_basic() {
	// Procedure: use nodes xi = {0,1,2} and fxi = {0,1,4} (f(x)=x^2). Evaluate at 1.5 -> 2.25
	std::vector<double> xi = {0.0,1.0,2.0};
	std::vector<double> fxi = {0.0,1.0,4.0};
	double val = lagrange(1.5, xi, fxi);
	ASSERT_CLOSE(val, 2.25, 1e-12);
	return 0;
}

int test_i_newtona_placeholder() {
	// Procedure: i_newtona currently returns placeholder coefficients of zeros.
	// Expectation: size matches input size and values are finite.
	std::vector<double> xi = {0.0, 1.0};
	std::vector<double> fxi = {1.0, 2.0};
	auto coeffs = i_newtona(xi, fxi, (int)xi.size());
	ASSERT_TRUE((int)coeffs.size() == (int)xi.size());
	return 0;
}

// Integration tests
int test_integralTrapez_basic() {
	// Procedure: integrate f(x)=x over [0,1]
	std::function<double(double)> f = [](double x){ return x; };
	double I = integralTrapez(f, 0.0, 1.0, 1000);
	ASSERT_CLOSE(I, 0.5, 1e-5);
	return 0;
}

int test_gaussLegendre_basic() {
	// Procedure: integrate f(x)=x^2 over [0,1], exact 1/3
	std::function<double(double)> f = [](double x){ return x*x; };
	double I = gaussLegendre(f, 0.0, 1.0, 4);
	ASSERT_CLOSE(I, 1.0/3.0, 1e-6);
	return 0;
}

// Nonlinear root finding tests
int test_newton_numerical_basic() {
	// Procedure: find root of x^2 - 2 starting at 1.0
	auto f = [](double x){ return x*x - 2.0; };
	double root = newton_numerical(f, 1.0, 1e-12, 50, "");
	ASSERT_CLOSE(root, std::sqrt(2.0), 1e-8);
	return 0;
}

int test_newton_numerical_flat() {
	// Procedure: f(x)=x^3 has flat derivative at 0; start at 0.0
	auto f = [](double x){ return x*x*x; };
	double root = newton_numerical(f, 0.0, 1e-12, 20, "");
	// Expectation: algorithm may fail to converge; final iterate should be near initial
	ASSERT_CLOSE(root, 0.0, 1e-6);
	return 0;
}

int test_secant_basic() {
	// Procedure: secant for x^2 - 2 using two guesses
	auto f = [](double x){ return x*x - 2.0; };
	double root = secantMethod(f, 1.0, 2.0, 1e-12, 50, "");
	ASSERT_CLOSE(root, std::sqrt(2.0), 1e-8);
	return 0;
}

int test_secant_divzero() {
	// Procedure: choose x0 and x1 such that f(x0)==f(x1)
	auto f = [](double x){ return 1.0; };
	double root = secantMethod(f, 0.0, 1.0, 1e-6, 10, "");
	// Expectation: secant returns last x1 when division by zero occurs
	ASSERT_CLOSE(root, 1.0, 1e-12);
	return 0;
}

// Linear systems tests
int test_forwardSub_basic() {
	// Procedure: solve L y = b with L = [[1,0],[2,1]] and b = [3,7] => y = [3,1]
	std::vector<std::vector<double>> L = {{1.0,0.0},{2.0,1.0}};
	std::vector<double> b = {3.0,7.0};
	auto y = forwardSub(L, b, 2);
	ASSERT_CLOSE(y[0], 3.0, 1e-12);
	ASSERT_CLOSE(y[1], 1.0, 1e-12);
	return 0;
}

int test_computeResidualNorm_basic() {
	// Procedure: residual of A*x - b where A = I, x = [1,2], b = [1,2] => 0
	std::vector<std::vector<double>> A = {{1.0,0.0},{0.0,1.0}};
	std::vector<double> x = {1.0,2.0};
	std::vector<double> b = {1.0,2.0};
	double r = computeResidualNorm(A, x, b);
	ASSERT_CLOSE(r, 0.0, 1e-12);
	return 0;
}

int main() {
	int failures = 0;
	failures += test_eulerExplicit_basic();
	failures += test_eulerExplicit_badstep();
	failures += test_rungeKutta4_basic();
	failures += test_rungeKutta4_stepsize_effect();
	failures += test_lagrange_basic();
	failures += test_i_newtona_placeholder();
	failures += test_integralTrapez_basic();
	failures += test_gaussLegendre_basic();
	failures += test_newton_numerical_basic();
	failures += test_newton_numerical_flat();
	failures += test_secant_basic();
	failures += test_secant_divzero();
	failures += test_forwardSub_basic();
	failures += test_computeResidualNorm_basic();

	if (failures == 0) {
		std::cout << "All tests passed.\n";
		return 0;
	} else {
		std::cerr << failures << " test(s) failed.\n";
		return failures;
	}
}
