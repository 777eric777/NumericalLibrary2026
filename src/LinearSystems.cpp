#include "../include/LinearSystems.h"
#include <stdexcept>
#include <cmath>
#include <iostream>

namespace NumericalLibrary {

std::vector<std::vector<double>> gauss(int N, std::vector<std::vector<double>>& a, std::vector<double>& b) {
	int h = 0;
	int k = 0;
	while (h < N && k < N) {
		for (int i = h + 1; i < N; i++) {
			double f = a[i][h] / a[k][h];
			a[i][k] = 0;
			b[i] = b[i] - b[h] * f;
			for (int j = k + 1; j < N; j++) {
				a[i][j] = a[i][j] - a[h][j] * f;
			}
		}
		h++;
		k++;
	}
	return a;
}

std::vector<double> xGauss(const std::vector<std::vector<double>>& a, const std::vector<double>& b, int n) {
	std::vector<double> x(n);
	for (int i = n - 1; i >= 0; i--) {
		x[i] = b[i];
		for (int j = i + 1; j < n; j++) {
			x[i] -= a[i][j] * x[j];
		}
		x[i] /= a[i][i];
	}
	return x;
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> LU_decomposition(const std::vector<std::vector<double>>& A, int N) {
	std::vector<std::vector<double>> L(N, std::vector<double>(N, 0.0));
	std::vector<std::vector<double>> U(N, std::vector<double>(N, 0.0));

	for (int i = 0; i < N; ++i) {
		for (int j = i; j < N; ++j) {
			double sum = 0.0;
			for (int k = 0; k < i; ++k) sum += L[i][k] * U[k][j];
			U[i][j] = A[i][j] - sum;
		}
		L[i][i] = 1.0;
		for (int k = i + 1; k < N; ++k) {
			double sum = 0.0;
			for (int j = 0; j < i; ++j) sum += L[k][j] * U[j][i];
			if (fabs(U[i][i]) < 1e-18) throw std::runtime_error("Zero or tiny pivot in LU_decomposition; use pivoting");
			L[k][i] = (A[k][i] - sum) / U[i][i];
		}
	}
	return {L, U};
}

std::vector<double> forwardSub(const std::vector<std::vector<double>>& L, const std::vector<double>& b, int n) {
	std::vector<double> y(n);
	for (int i = 0; i < n; ++i) {
		double sum = 0.0;
		for (int j = 0; j < i; ++j) sum += L[i][j] * y[j];
		if (L[i][i] == 0.0) throw std::runtime_error("Zero pivot in L");
		y[i] = (b[i] - sum) / L[i][i];
	}
	return y;
}

double computeResidualNorm(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const std::vector<double>& b) {
	int n = (int)A.size();
	double maxErr = 0.0;
	for (int i = 0; i < n; ++i) {
		double s = 0.0;
		for (int j = 0; j < n; ++j) s += A[i][j] * x[j];
		double err = fabs(s - b[i]);
		if (err > maxErr) maxErr = err;
	}
	return maxErr;
}

} // namespace NumericalLibrary
