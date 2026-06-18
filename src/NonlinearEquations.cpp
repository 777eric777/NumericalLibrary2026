#include "../include/NonlinearEquations.h"
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>

namespace fs = std::filesystem;

namespace NumericalLibrary {

// Per-iteration record used for CSV diagnostics.
struct IterData { std::string method; double start0, start1; int iter; double x; double fx; double dx; };

// Append rows to a CSV file; create parent directories and header when needed.
static void writeIterationsToCSV(const std::string& filename, const std::vector<IterData>& rows) {
	if (filename.empty()) return; // nothing to do when caller didn't request output
	fs::path p(filename);
	if (p.has_parent_path()) fs::create_directories(p.parent_path());
	bool exists = fs::exists(p);
	std::ofstream f(filename, std::ios::app);
	if (!exists) f << "method,start0,start1,iter,x,f_x,dx\n";
	for (auto &r : rows) {
		f << r.method << "," << r.start0 << "," << r.start1 << "," << r.iter << "," << std::setprecision(15) << r.x << "," << r.fx << "," << r.dx << "\n";
	}
}

// Append a short run summary to a CSV file for result aggregation.
static void appendSummaryCSV(const std::string& filename, const std::string& funcName, const std::string& method, double start0, double start1, double root, double froot, int iters, bool converged) {
	if (filename.empty()) return;
	fs::path p(filename);
	if (p.has_parent_path()) fs::create_directories(p.parent_path());
	bool exists = fs::exists(p);
	std::ofstream f(filename, std::ios::app);
	if (!exists) f << "function,method,start0,start1,root,f_root,iters,converged\n";
	f << funcName << "," << method << "," << start0 << "," << start1 << "," << std::setprecision(15) << root << "," << froot << "," << iters << "," << (converged ? "1" : "0") << "\n";
}

/**
 * Numerical derivative (central difference).
 * Approximates f'(x) using (f(x+h)-f(x-h)) / (2h). Default h=1e-6.
 */
double derivativeNumerical(const std::function<double(double)>& f, double x, double h) {
	return (f(x + h) - f(x - h)) / (2.0 * h);
}

/**
 * Newton's method with analytical derivative.
 */
double newton_analytical(std::function<double(double)> df_analytical, std::function<double(double)> f, double x0, double tol, int max_iter, const std::string& outFile) {
	double x = x0;
	std::vector<IterData> rows;
	for (int i = 0; i < max_iter; ++i) {
		double fx = f(x);
		double dfx = df_analytical(x);
		if (std::fabs(dfx) < 1e-18) {
			std::cerr << "Derivative too small; aborting Newton (analytic)." << std::endl;
			break;
		}
		double x_next = x - fx / dfx;
		rows.push_back({ "newton_analytic", x0, NAN, i+1, x_next, f(x_next), std::fabs(x_next - x) });
		if (std::fabs(x_next - x) < tol) {
			writeIterationsToCSV(outFile, rows);
			appendSummaryCSV("data/summary.csv", "f", "newton_analytic", x0, NAN, x_next, f(x_next), i+1, true);
			return x_next;
		}
		x = x_next;
	}
	writeIterationsToCSV(outFile, rows);
	appendSummaryCSV("data/summary.csv", "f", "newton_analytic", x0, NAN, x, f(x), max_iter, false);
	return x;
}

/**
 * Newton's method with numerical derivative (central difference).
 */
double newton_numerical(std::function<double(double)> f, double x0, double tol, int max_iter, const std::string& outFile) {
	double x = x0;
	std::vector<IterData> rows;
	for (int i = 0; i < max_iter; ++i) {
		double fx = f(x);
		double dfx = derivativeNumerical(f, x, 1e-6);
		if (std::fabs(dfx) < 1e-18) {
			std::cerr << "Derivative too small; aborting Newton (numeric)." << std::endl;
			break;
		}
		double x_next = x - fx / dfx;
		rows.push_back({ "newton_numeric", x0, NAN, i+1, x_next, f(x_next), std::fabs(x_next - x) });
		if (std::fabs(x_next - x) < tol) {
			writeIterationsToCSV(outFile, rows);
			appendSummaryCSV("data/summary.csv", "f", "newton_numeric", x0, NAN, x_next, f(x_next), i+1, true);
			return x_next;
		}
		x = x_next;
	}
	writeIterationsToCSV(outFile, rows);
	appendSummaryCSV("data/summary.csv", "f", "newton_numeric", x0, NAN, x, f(x), max_iter, false);
	return x;
}

/**
 * Secant method.
 */
double secantMethod(std::function<double(double)> f, double x0, double x1, double tolerance, int maxIterations, const std::string& outFile) {
	double x2 = x1;
	std::vector<IterData> rows;
	for (int i = 0; i < maxIterations; ++i) {
		double fx0 = f(x0);
		double fx1 = f(x1);
		if (std::abs(fx1 - fx0) < 1e-18) {
			writeIterationsToCSV(outFile, rows);
			appendSummaryCSV("data/summary.csv", "f", "secant", x0, x1, x1, f(x1), i, false);
			std::cout << "Division by zero in secant; returning last x1." << std::endl;
			return x1;
		}
		x2 = x1 - (fx1 * (x1 - x0)) / (fx1 - fx0);
		rows.push_back({ "secant", x0, x1, i+1, x2, f(x2), std::fabs(x2 - x1) });
		if (std::abs(f(x2)) < tolerance) {
			writeIterationsToCSV(outFile, rows);
			appendSummaryCSV("data/summary.csv", "f", "secant", x0, x1, x2, f(x2), i+1, true);
			return x2;
		}
		x0 = x1;
		x1 = x2;
	}
	writeIterationsToCSV(outFile, rows);
	appendSummaryCSV("data/summary.csv", "f", "secant", x0, x1, x2, f(x2), maxIterations, false);
	std::cout << "Secant: max iterations reached; returning last approximation." << std::endl;
	return x2;
}

/**
 * Regula falsi (false position) method.
 */
bool false_position(const std::function<double(double)>& f, double a, double b, double& root,
	double tol, int maxIter, bool verbose,
	std::vector<std::pair<int,double>>* iter_log)
{
	if (a > b) std::swap(a, b);
	double fa = f(a);
	double fb = f(b);
	if (std::fabs(fa) < tol) { root = a; return true; }
	if (std::fabs(fb) < tol) { root = b; return true; }
	if (fa * fb > 0) {
		std::cout << "No sign change in [" << a << "," << b << "] for false_position." << std::endl;
		return false;
	}
	double c = a;
	double prev_c = c;
	for (int iter = 0; iter < maxIter; ++iter) {
		c = (a * fb - b * fa) / (fb - fa);
		double fc = f(c);
		if (iter_log) iter_log->emplace_back(iter, c);
		if (verbose) {
			std::cout << "falsi iter=" << iter << " a=" << a << " b=" << b << " c=" << c << " f(c)=" << fc << "\n";
		}
		prev_c = c;
		if (std::fabs(fc) < tol || std::fabs(b - a) < tol) {
			root = c;
			return true;
		}
		if (fa * fc < 0) {
			b = c; fb = fc;
		} else {
			a = c; fa = fc;
		}
	}
	root = c;
	std::cout << "false_position: max iterations reached; returning approximation " << std::setprecision(15) << root << std::endl;
	return true;
}

/**
 * Scan [a,b] by subintervals and apply false_position where sign changes
 * occur. Returns distinct roots and optional per-root logs.
 */
bool find_roots_false_position(const std::function<double(double)>& f, double a, double b, int l_podpodzialow,
	std::vector<double>& roots, bool verbose, double tol, int maxIter,
	std::vector<std::vector<std::pair<int,double>>>* all_iter_logs)
{
	roots.clear();
	if (a > b) std::swap(a, b);
	if (l_podpodzialow < 1) l_podpodzialow = 1;
	double step = (b - a) / l_podpodzialow;
	std::vector<std::pair<double, std::vector<std::pair<int,double>>>> found;
	for (int i = 0; i < l_podpodzialow; ++i) {
		double left = a + i * step;
		double right = (i == l_podpodzialow - 1) ? b : left + step;
		double fl = f(left);
		double fr = f(right);
		if (std::fabs(fl) < tol) found.emplace_back(left, std::vector<std::pair<int,double>>());
		if (std::fabs(fr) < tol) found.emplace_back(right, std::vector<std::pair<int,double>>());
		double mid = 0.5 * (left + right);
		double fm = f(mid);
		if (std::fabs(fm) < tol) found.emplace_back(mid, std::vector<std::pair<int,double>>());
		if (fl * fr <= 0 && !(std::fabs(fl) < tol || std::fabs(fr) < tol)) {
			double r;
			std::vector<std::pair<int,double>> log;
			if (false_position(f, left, right, r, tol, maxIter, verbose, &log)) {
				found.emplace_back(r, std::move(log));
			}
		}
	}
	if (found.empty()) return false;
	const double mergeTol = std::max(tol, 1e-10);
	std::sort(found.begin(), found.end(), [](auto &A, auto &B){ return A.first < B.first; });
	std::vector<double> uniq;
	std::vector<std::vector<std::pair<int,double>>> uniqLogs;
	for (auto &pr : found) {
		if (uniq.empty() || std::fabs(pr.first - uniq.back()) > mergeTol) {
			uniq.push_back(pr.first);
			uniqLogs.push_back(pr.second);
		}
	}
	roots.swap(uniq);
	if (all_iter_logs) *all_iter_logs = std::move(uniqLogs);
	return !roots.empty();
}

/**
 * Classic bisection method on [a,b].
 */
bool bisection(const std::function<double(double)>& f, double a, double b, double& root,
	double tol, int maxIter, bool verbose,
	std::vector<std::pair<int,double>>* iter_log)
{
	if (a > b) std::swap(a, b);
	double fa = f(a);
	double fb = f(b);
	if (std::fabs(fa) < tol) { root = a; return true; }
	if (std::fabs(fb) < tol) { root = b; return true; }
	if (fa * fb > 0) {
		std::cout << "bisection: no sign change in interval." << std::endl;
		return false;
	}
	double c = a;
	double prev_c = c;
	for (int iter = 0; iter < maxIter; ++iter) {
		c = 0.5 * (a + b);
		double fc = f(c);
		if (iter_log) iter_log->emplace_back(iter, c);
		if (verbose) {
			std::cout << "bisection iter=" << iter << " a=" << a << " b=" << b << " c=" << c << " f(c)=" << fc << "\n";
		}
		if (std::fabs(fc) < tol || 0.5 * (b - a) < tol) { root = c; return true; }
		if (fa * fc < 0) { b = c; fb = fc; } else { a = c; fa = fc; }
		prev_c = c;
	}
	root = c;
	return true;
}

/**
 * Find roots using bisection across subdivided intervals.
 */
bool find_roots_bisection(const std::function<double(double)>& f, double a, double b, int l_podpodzialow,
	std::vector<double>& roots, bool verbose, double tol, int maxIter,
	std::vector<std::vector<std::pair<int,double>>>* all_iter_logs)
{
	roots.clear();
	if (a > b) std::swap(a, b);
	if (l_podpodzialow < 1) l_podpodzialow = 1;
	double step = (b - a) / l_podpodzialow;
	std::vector<std::pair<double, std::vector<std::pair<int,double>>>> found;
	for (int i = 0; i < l_podpodzialow; ++i) {
		double left = a + i * step;
		double right = (i == l_podpodzialow - 1) ? b : left + step;
		double fl = f(left);
		double fr = f(right);
		if (std::fabs(fl) < tol) found.emplace_back(left, std::vector<std::pair<int,double>>());
		if (std::fabs(fr) < tol) found.emplace_back(right, std::vector<std::pair<int,double>>());
		double mid = 0.5 * (left + right);
		double fm = f(mid);
		if (std::fabs(fm) < tol) found.emplace_back(mid, std::vector<std::pair<int,double>>());
		if (fl * fr <= 0 && !(std::fabs(fl) < tol || std::fabs(fr) < tol)) {
			double r;
			std::vector<std::pair<int,double>> log;
			if (bisection(f, left, right, r, tol, maxIter, verbose, &log)) found.emplace_back(r, std::move(log));
		}
	}
	if (found.empty()) return false;
	const double mergeTol = std::max(tol, 1e-10);
	std::sort(found.begin(), found.end(), [](auto &A, auto &B){ return A.first < B.first; });
	std::vector<double> uniq;
	std::vector<std::vector<std::pair<int,double>>> uniqLogs;
	for (auto &pr : found) {
		if (uniq.empty() || std::fabs(pr.first - uniq.back()) > mergeTol) {
			uniq.push_back(pr.first);
			uniqLogs.push_back(pr.second);
		}
	}
	roots.swap(uniq);
	if (all_iter_logs) *all_iter_logs = std::move(uniqLogs);
	return !roots.empty();
}

} // namespace NumericalLibrary
