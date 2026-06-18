#include "NumericalLibrary.h"
#include <iostream>

int main() {
	using namespace NumericalLibrary;

	// Example: Newton analytic for f(x) = x^2 - 2
	auto f = [](double x) { return x * x - 2.0; };
	auto df = [](double x) { return 2.0 * x; };
	double root = newton_analytical(df, f, 1.0, 1e-12, 50, "");
	std::cout << "root = " << root << "\n";

	// Example: RK4 for y' = -y, y(0)=1
	auto rhs = [](double t, double y) { return -y; };
	auto [t, y] = rungeKutta4(1.0, 0.0, 2.0, 0.1, rhs);
	std::cout << "RK4 last value y(" << t.back() << ") = " << y.back() << "\n";

	return 0;
}
