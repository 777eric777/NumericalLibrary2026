#include "../include/Interpolation.h"
#include <cmath>

namespace NumericalLibrary {

double lagrange(double x, const std::vector<double>& xi, const std::vector<double>& fxi) {
	double wynik = 0.0;
	for (size_t i = 0; i < xi.size(); ++i) {
		double li = 1.0;
		for (size_t j = 0; j < xi.size(); ++j) {
			if (j != i) li *= (x - xi[j]) / (xi[i] - xi[j]);
		}
		wynik += fxi[i] * li;
	}
	return wynik;
}

std::vector<double> i_newtona(const std::vector<double>& xi, const std::vector<double>& fxi, int size) {
	std::vector<double> wsp;
	// Placeholder: actual implementation of divided differences
	wsp.reserve(size);
	for (int i = 0; i < size; ++i) wsp.push_back(0.0);
	return wsp;
}

} // namespace NumericalLibrary
