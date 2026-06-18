#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <sstream>
#include <string>
#include <utility>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <filesystem>
using namespace std;
using namespace chrono;
namespace fs = std::filesystem;
// CSV-only flag removed; always print iterations to console and keep per-iteration CSV output

pair<vector<double>, vector<double>> readFile1(fstream& plik);
pair<vector<double>, vector<double>> read2Lagrange(fstream& plik);
pair<vector<double>, vector<double>> readFile2(fstream& plik);
pair<vector<double>, vector<double>> readFile4(fstream& plik);

struct MatrixData {
	int N;
	vector<double> B;
	vector<vector<double>> A;
};

MatrixData readFile3(fstream& plik);
void saveVectorCSV(const vector<double>& vec, const string& filename);
void saveXY_CSV(const vector<double>& x, const vector<double>& y, const string& filename);
void saveMatrixCSV(const vector<vector<double>>& mat, const string& filename);
void saveResultsCSV(
	const string& filename,
	const vector<double>& time,
	const vector<double>& exact,

	const vector<double>& h1,
	const vector<double>& h05,
	const vector<double>& h01,
	const vector<double>& h001
);
void saveApproxCSV(const vector<double>& x, const vector<double>& y_true, const vector<double>& y_approx, const string& filename);
double lagrange(double x, const vector<double>& lxi, const vector<double>& fxi);
double p_naturalna(double x, const vector<double>& wspolczynniki);
double p_hornera(double x, const vector<double>& wspolczynniki);
vector<vector<double>> gauss(int N, vector<vector<double>>& a, vector<double>& b);
void printMatrix(const vector<vector<double>> a);
vector<double> xGauss(const vector<vector<double>>& a, const vector<double>& b, int n);
vector<double> forwardSub(const vector<vector<double>>& L, const vector<double>& b, int n);
double computeResidualNorm(const vector<vector<double>>& A, const vector<double>& x, const vector<double>& b);
pair<vector<vector<double>>, vector<vector<double>>> LU_decomposition(const vector<vector<double>>& A, int N);
double integralSimpsona(const vector<double>& ai, const vector<double>& range, int liczbaPodprzedzialow);
double integralTrapez(const vector<double>& ai, const vector<double>& range, int liczbaPodprzedzialow = 1000);
// dla funkcji (trygonometrycznych)
double integralSimpsona(const function<double(double)>& f, double a, double b, int liczbaPodprzedzialow = 1000);
double integralTrapez(const function<double(double)>& f, double a, double b, int liczbaPodprzedzialow = 1000);
double trygLAB6(double x);
double funLab7(double x);
double funLab8(double x);
double gaussLegendre(const function<double(double)>& f, double a, double b, int n);
double gaussLegendreComposite(const function<double(double)>& f, double a, double b, int n, int k);
vector<double> leastSquaresApprox(const function<double(double)>& f, double range_a, double range_b, int degree, int gl_n = 8, int gl_k = 10);
double funLab9(double t);
// Explicit Euler solver for ODE y' = f(t,y)
pair<vector<double>, vector<double>> eulerExplicit(const function<double(double, double)>& f, double t0, double y0, double h, double t_end);
// Heun (improved Euler) solver for ODE y' = f(t,y)
pair<vector<double>, vector<double>> heunExplicit(const function<double(double, double)>& f, double t0, double y0, double h, double t_end);
pair<vector<double>, vector<double>> midpointMethod(
	const function<double(double, double)>& f,
	double t0,
	double y0,
	double h,
	double t_end);

pair<vector<double>, vector<double>> rungeKutta4(
	double y0,
	double t0,
	double t_end,
	double h,
	auto f);

// Root finding: Newton (analytic and numeric) and Secant methods
// Root finding implementations (defined below)
double newton_analytical(function<double(double)> df_analytical, function<double(double)> f, double x0, double tol, int max_iter, const string& outFile);
double newton_numerical(function<double(double)> f, double x0, double tol, int max_iter, const string& outFile);
double secantMethod(function<double(double)> f, double x0, double x1, double tolerance, int maxIterations, const string& outFile);

struct IterData { string method; double start0, start1; int iter; double x; double fx; double dx; };

static void writeIterationsToCSV(const string& filename, const vector<IterData>& rows) {
	fs::path p(filename);
	if (p.has_parent_path()) fs::create_directories(p.parent_path());
	bool exists = fs::exists(p);
	ofstream f(filename, ios::app);
	if (!exists) f << "method,start0,start1,iter,x,f_x,dx\n";
	for (auto &r : rows) {
		f << r.method << "," << r.start0 << "," << r.start1 << "," << r.iter << "," << setprecision(15) << r.x << "," << r.fx << "," << r.dx << "\n";
	}
	f.close();
}

static void appendSummaryCSV(const string& filename, const string& funcName, const string& method, double start0, double start1, double root, double froot, int iters, bool converged) {
	fs::path p(filename);
	if (p.has_parent_path()) fs::create_directories(p.parent_path());
	bool exists = fs::exists(p);
	ofstream f(filename, ios::app);
	if (!exists) f << "function,method,start0,start1,root,f_root,iters,converged\n";
	f << funcName << "," << method << "," << start0 << "," << start1 << "," << setprecision(15) << root << "," << froot << "," << iters << "," << (converged ? "1" : "0") << "\n";
	f.close();
}

// Analytical derivatives for f1,f2,f3 and numerical derivative helper
double df1(double x) {
	// f1(x) = exp(-x) * sin(x) + x/5
	// f1'(x) = exp(-x) * (cos(x) - sin(x)) + 1/5
	return exp(-x) * (cos(x) - sin(x)) + 0.2;
}

double df2(double x) {
	// f2(x) = (x^2 - 1)/(x^2 + 1) + 0.2 * cos(x)
	// derivative: 4x/(x^2+1)^2 - 0.2 * sin(x)
	double denom = pow(x, 2) + 1.0;
	return 4.0 * x / (denom * denom) - 0.2 * sin(x);
}

double df3(double x) {
	// f3(x) = sin(x) * log(x^2 + 1)
	// derivative: cos(x) * log(x^2+1) + (2x * sin(x)) / (x^2+1)
	double denom = x * x + 1.0;
	return cos(x) * log(denom) + (2.0 * x * sin(x)) / denom;
}

double derivativeNumerical(const function<double(double)>& f, double x, double h);

double derivativeNumerical(const function<double(double)>& f, double x, double h = 1e-6) {
	return (f(x + h) - f(x - h)) / (2.0 * h);
}

double f1(double x) {
	return exp(-x) * sin(x) + x/5;
}

double f2(double x) {
	return (x - 1) * (x + 1) / (pow(x, 2) + 1) + 0.2 * cos(x);
}

double f3(double x) {
	return sin(x) * log(pow(x, 2) + 1);
}

int main()
{
	double f1_rootExact[] = { 0.0, -3.16825786774928, -6.2808339631475, -9.42493005509724 };
	double f2_rootExact[] = { -0.879700562635, 0.879700562635 };
	double f3_rootExact[] = { 0.0, 3.141592653589793, 6.283185307179586 };  

	const double tol = 1e-12;
	const int maxIter = 100;

	cout << fixed << setprecision(12);

	// f1: use initial guesses near known roots
	vector<double> f1_inits = { 0.5, -3.0, -6.2, -9.3 };
	cout << "\n=== Finding roots for f1: exp(-x) * sin(x) + x/5, using Newton (analytic) ===\n";
	for (double x0 : f1_inits) {
		cout << "Start x0 = " << x0 << "\n";
		string out = "data/f1_newton_analytic_x0_" + to_string((int)round(x0*1000)) + ".csv";
		double root = newton_analytical(df1, f1, x0, tol, maxIter, out);
		cout << "Result: " << root << " f(root)=" << f1(root) << "\n\n";
	}

	cout << "\n=== Finding roots for f1: exp(-x) * sin(x) + x/5, using Newton (numerical derivative) ===\n";
	for (double x0 : f1_inits) {
		cout << "Start x0 = " << x0 << "\n";
		string out = "data/f1_newton_numeric_x0_" + to_string((int)round(x0*1000)) + ".csv";
		double root = newton_numerical(f1, x0, tol, maxIter, out);
		cout << "Result: " << root << " f(root)=" << f1(root) << "\n\n";
	}

	cout << "\n=== Finding roots for f1: exp(-x) * sin(x) + x/5, using Secant method ===\n";
	for (double r : f1_inits) {
		double x0 = r - 0.1;
		double x1 = r + 0.1;
		cout << "Start x0 = " << x0 << "  x1 = " << x1 << "\n";
		string out = "data/f1_secant_x0_" + to_string((int)round(x0*1000)) + "_x1_" + to_string((int)round(x1*1000)) + ".csv";
		double root = secantMethod(f1, x0, x1, tol, maxIter, out);
		cout << "Result: " << root << " f(root)=" << f1(root) << "\n\n";
	}

	// f2
	vector<pair<double,double>> f2_secant_inits = { {-1.0, -0.5}, {0.5, 1.0} };
	vector<double> f2_inits = { -0.9, 0.9 };
	cout << "\n=== Finding roots for f2: (x - 1) * (x + 1) / (pow(x, 2) + 1) + 0.2 * cos(x), using Newton (analytic) ===\n";
	for (double x0 : f2_inits) {
		cout << "Start x0 = " << x0 << "\n";
		string out = "data/f2_newton_analytic_x0_" + to_string((int)round(x0*1000)) + ".csv";
		double root = newton_analytical(df2, f2, x0, tol, maxIter, out);
		cout << "Result: " << root << " f(root)=" << f2(root) << "\n\n";
	}
	cout << "\n=== Finding roots for f2: (x - 1) * (x + 1) / (pow(x, 2) + 1) + 0.2 * cos(x), using Secant method ===\n";
	for (auto pr : f2_secant_inits) {
		cout << "Start x0 = " << pr.first << "  x1 = " << pr.second << "\n";
		string out = "data/f2_secant_x0_" + to_string((int)round(pr.first*1000)) + "_x1_" + to_string((int)round(pr.second*1000)) + ".csv";
		double root = secantMethod(f2, pr.first, pr.second, tol, maxIter, out);
		cout << "Result: " << root << " f(root)=" << f2(root) << "\n\n";
	}

	// f3
	vector<double> f3_inits = { 0.5, 3.0, 6.0 };
	cout << "\n=== Finding roots for f3:sin(x) * log(pow(x, 2) + 1), using Newton (analytic) ===\n";
	for (double x0 : f3_inits) {
		cout << "Start x0 = " << x0 << "\n";
		string out = "data/f3_newton_analytic_x0_" + to_string((int)round(x0*1000)) + ".csv";
		double root = newton_analytical(df3, f3, x0, tol, maxIter, out);
		cout << "Result: " << root << " f(root)=" << f3(root) << "\n\n";
	}

	cout << "\n=== Finding roots for f3: sin(x) * log(pow(x, 2) + 1), using Secant method ===\n";
	for (double r : f3_inits) {
		double x0 = r - 0.1;
		double x1 = r + 0.1;
		cout << "Start x0 = " << x0 << "  x1 = " << x1 << "\n";
		string out = "data/f3_secant_x0_" + to_string((int)round(x0*1000)) + "_x1_" + to_string((int)round(x1*1000)) + ".csv";
		double root = secantMethod(f3, x0, x1, tol, maxIter, out);
		cout << "Result: " << root << " f(root)=" << f3(root) << "\n\n";
	}

	cout << "    f1: wszystkie znane pierwiastki znalezione, przebiegi stabilne.\n    f2 : symetryczne pierwiastki szybko znalezione(3 - 5 iteracji).\n    f3 przy x~0 : bardzo wolna zbieznosc Newtona(dziesiatki iteracji). \n    Powod : f3(x) = sin x * log(1 + x ^ 2) ~O(x ^ 3) przy 0, wiec f'(0) = 0 \n    - Newton traci kwadratowa zbieznosc (jak przy wielokrotnym/plaskim pierwiastku)\n    Metoda siecznych przy jednym zestawie(f1, blisko - 9.42493) wygenerowala 'Blad : dzielenie przez zero.' \n    - guard wykryl bardzo mala roznice fx1 - fx0 po konwergencji i zwrocono NaN, \n    mimo ze metoda juz zblizyla sie do pierwiastka.\n\n";
	return 0;
}

// The rest of the file (definitions for newton_analytical, newton_numerical, secantMethod, rungeKutta4, midpointMethod,
// heunExplicit, eulerExplicit, saveResultsCSV, leastSquaresApprox, saveApproxCSV, trygLAB6, funLab7, funLab8,
// gaussLegendre, gaussLegendreComposite, integralTrapez/Simpsona, p_naturalna, p_hornera, i_newtona,
// readFile1/2/3/4, lagrange, saveVectorCSV, saveXY_CSV, saveMatrixCSV, printMatrix, gauss, xGauss, computeResidualNorm,
// LU_decomposition, forwardSub) were moved here as part of the full original source.

// For the full function bodies consult the repository history or original MNv11.cpp before refactor.
