#include "../include/ODE.h"
#include <cmath>

namespace NumericalLibrary {

std::pair<std::vector<double>, std::vector<double>> rungeKutta4(double y0, double t0, double t_end, double h, const std::function<double(double,double)>& f) {
	std::vector<double> time;
	std::vector<double> values;

	double t = t0;
	double y = y0;

	time.push_back(t);
	values.push_back(y);

	while (t < t_end) {
		double k1 = f(t, y);
		double k2 = f(t + h / 2.0, y + (h / 2.0) * k1);
		double k3 = f(t + h / 2.0, y + (h / 2.0) * k2);
		double k4 = f(t + h, y + h * k3);
		y = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		t = t + h;
		time.push_back(t);
		values.push_back(y);
	}

	return {time, values};
}

std::pair<std::vector<double>, std::vector<double>> midpointMethod(const std::function<double(double, double)>& f, double t0, double y0, double h, double t_end) {
	std::vector<double> time;
	std::vector<double> temperature;

	double t = t0;
	double y = y0;

	time.push_back(t);
	temperature.push_back(y);

	while (t < t_end) {
		double k1 = f(t, y);
		double ymid = y + (h / 2.0) * k1;
		double tmid = t + h / 2.0;
		double k2 = f(tmid, ymid);
		y = y + h * k2;
		t = t + h;
		time.push_back(t);
		temperature.push_back(y);
	}

	return {time, temperature};
}

std::pair<std::vector<double>, std::vector<double>> heunExplicit(const std::function<double(double, double)>& f, double t0, double y0, double h, double t_end) {
	std::vector<double> t_values;
	std::vector<double> y_values;

	double t = t0;
	double y = y0;

	int steps = static_cast<int>((t_end - t0) / h);

	for (int i = 0; i <= steps; ++i) {
		t_values.push_back(t);
		y_values.push_back(y);
		double k1 = f(t, y);
		double y_predict = y + h * k1;
		double k2 = f(t + h, y_predict);
		y += (h / 2.0) * (k1 + k2);
		t += h;
	}

	return {t_values, y_values};
}

std::pair<std::vector<double>, std::vector<double>> eulerExplicit(const std::function<double(double, double)>& f, double t0, double y0, double h, double t_end) {
	std::vector<double> t_values;
	std::vector<double> y_values;

	double t = t0;
	double y = y0;

	int steps = static_cast<int>((t_end - t0) / h);

	for (int i = 0; i <= steps; ++i) {
		t_values.push_back(t);
		y_values.push_back(y);
		double dy_dt = f(t, y);
		y += h * dy_dt;
		t += h;
	}

	return {t_values, y_values};
}

} // namespace NumericalLibrary
