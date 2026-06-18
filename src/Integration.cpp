#include "../include/Integration.h"
#include <cmath>

namespace NumericalLibrary {

double gaussLegendre(const std::function<double(double)>& f, double a, double b, int n) {
	if (n < 2) n = 2;
	if (n > 8) n = 8;

	static const double t2[2] = { -0.5773502691896257, 0.5773502691896257 };
	static const double w2[2] = { 1.0, 1.0 };
	static const double t3[3] = { -0.7745966692414834, 0.0, 0.7745966692414834 };
	static const double w3[3] = { 0.5555555555555556, 0.8888888888888888, 0.5555555555555556 };
	static const double t4[4] = { -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526 };
	static const double w4[4] = { 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538 };
	static const double t5[5] = { -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640 };
	static const double w5[5] = { 0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891 };
	static const double t6[6] = { -0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.9324695142031521 };
	static const double w6[6] = { 0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704 };
	static const double t7[7] = { -0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585 };
	static const double w7[7] = { 0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.2797053914892766, 0.1294849661688697 };
	static const double t8[8] = {
		-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
		 0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363
	};
	static const double w8[8] = {
		0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
		0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763
	};

	double c1 = 0.5 * (b - a);
	double c2 = 0.5 * (b + a);
	double sum = 0.0;

	if (n == 2) {
		for (int i = 0; i < 2; ++i) {
			double x = c1 * t2[i] + c2;
			sum += w2[i] * f(x);
		}
	}
	else if (n == 3) {
		for (int i = 0; i < 3; ++i) {
			double x = c1 * t3[i] + c2;
			sum += w3[i] * f(x);
		}
	}
	else if (n == 4) {
		for (int i = 0; i < 4; ++i) {
			double x = c1 * t4[i] + c2;
			sum += w4[i] * f(x);
		}
	}
	else if (n == 5) {
		for (int i = 0; i < 5; ++i) {
			double x = c1 * t5[i] + c2;
			sum += w5[i] * f(x);
		}
	}
	else if (n == 6) {
		for (int i = 0; i < 6; ++i) {
			double x = c1 * t6[i] + c2;
			sum += w6[i] * f(x);
		}
	}
	else if (n == 7) {
		for (int i = 0; i < 7; ++i) {
			double x = c1 * t7[i] + c2;
			sum += w7[i] * f(x);
		}
	}
	else if (n == 8) {
		for (int i = 0; i < 8; ++i) {
			double x = c1 * t8[i] + c2;
			sum += w8[i] * f(x);
		}
	}

	return c1 * sum;
}


double gaussLegendreComposite(const std::function<double(double)>& f, double a, double b, int n, int k) {
	if (a == b) return 0.0;
	if (a > b) std::swap(a, b);
	if (k < 1) k = 1;
	double h = (b - a) / k;
	double total = 0.0;
	for (int i = 0; i < k; ++i) {
		double sub_a = a + i * h;
		double sub_b = sub_a + h;
		total += gaussLegendre(f, sub_a, sub_b, n);
	}
	return total;
}


double integralTrapez(const std::vector<double>& ai, const std::vector<double>& range, int liczbaPodprzedzialow) {
	if (ai.empty() || range.size() < 2) return 0.0;
	double a = range[0];
	double b = range[1];
	if (a == b) return 0.0;
	if (a > b) std::swap(a, b);
	int n = std::max(1, liczbaPodprzedzialow);
	double h = (b - a) / n;
	double s = 0.5 * (/*p_hornera(a, ai)*/ 0.0 + /*p_hornera(b, ai)*/ 0.0);
	for (int i = 1; i < n; ++i) s += /*p_hornera(a + i * h, ai)*/ 0.0;
	return s * h;
}


double integralSimpsona(const std::vector<double>& ai, const std::vector<double>& range, int liczbaPodprzedzialow) {
	if (ai.empty() || range.size() < 2) return 0.0;
	double a = range[0];
	double b = range[1];
	if (a == b) return 0.0;
	if (a > b) std::swap(a, b);
	int m = liczbaPodprzedzialow;
	if (m % 2 == 1) --m;
	double h = (b - a) / m;
	double sum = 0.0;
	sum += /*p_hornera(a, ai)*/ 0.0;
	sum += /*p_hornera(b, ai)*/ 0.0;
	for (int i = 1; i < m; ++i) {
		double x = a + i * h;
		double fx = /*p_hornera(x, ai)*/ 0.0;
		if (i % 2 == 0) sum += 2.0 * fx;
		else sum += 4.0 * fx;
	}
	return (h / 3.0) * sum;
}


double integralTrapez(const std::function<double(double)>& f, double a, double b, int liczbaPodprzedzialow) {
	if (a == b) return 0.0;
	if (a > b) std::swap(a, b);
	int n = std::max(1, liczbaPodprzedzialow);
	double h = (b - a) / n;
	double s = 0.5 * (f(a) + f(b));
	for (int i = 1; i < n; ++i) s += f(a + i * h);
	return s * h;
}


double integralSimpsona(const std::function<double(double)>& f, double a, double b, int liczbaPodprzedzialow) {
	if (a == b) return 0.0;
	if (a > b) std::swap(a, b);
	int n = std::max(2, liczbaPodprzedzialow);
	if (n % 2 == 1) --n;
	double h = (b - a) / n;
	double sum = f(a) + f(b);
	for (int i = 1; i < n; ++i) {
		double x = a + i * h;
		if (i % 2 == 0) sum += 2.0 * f(x);
		else sum += 4.0 * f(x);
	}
	return (h / 3.0) * sum;
}

} // namespace NumericalLibrary
