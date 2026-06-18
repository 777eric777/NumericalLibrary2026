#include "../include/Approximation.h"
#include "../include/Integration.h"
#include <stdexcept>
#include <cmath>

namespace NumericalLibrary {

std::vector<double> leastSquaresApprox(const std::function<double(double)>& f, double range_a, double range_b, int degree, int gl_n, int gl_k) {
	int m = degree + 1;
	std::vector<std::vector<double>> M(m, std::vector<double>(m, 0.0));
	std::vector<double> rhs(m, 0.0);

	auto integrand_M = [&](int i, int j, double x) {
		return pow(x, i + j);
	};
	auto integrand_rhs = [&](int i, double x) {
		return f(x) * pow(x, i);
	};

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			auto g = [&](double x)->double { return integrand_M(i, j, x); };
			M[i][j] = gaussLegendreComposite(g, range_a, range_b, gl_n, gl_k);
		}
		auto h = [&](double x)->double { return integrand_rhs(i, x); };
		rhs[i] = gaussLegendreComposite(h, range_a, range_b, gl_n, gl_k);
	}

	int N = m;
	for (int k = 0; k < N; ++k) {
		int piv = k;
		for (int i = k + 1; i < N; ++i) if (fabs(M[i][k]) > fabs(M[piv][k])) piv = i;
		if (piv != k) { std::swap(M[piv], M[k]); std::swap(rhs[piv], rhs[k]); }
		double pivot = M[k][k];
		if (fabs(pivot) < 1e-18) throw std::runtime_error("Macierz normalna jest osobliwa w leastSquaresApprox");
		for (int j = k; j < N; ++j) M[k][j] /= pivot;
		rhs[k] /= pivot;
		for (int i = 0; i < N; ++i) if (i != k) {
			double fac = M[i][k];
			for (int j = k; j < N; ++j) M[i][j] -= fac * M[k][j];
			rhs[i] -= fac * rhs[k];
		}
	}

	return rhs;
}

} // namespace NumericalLibrary
