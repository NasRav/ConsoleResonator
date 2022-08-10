#pragma once

#include <complex>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "Resonator.h"

using namespace std;

class Linear_2D : public Resonator
{
public:
	Linear_2D();
	Linear_2D(int, int);
	Linear_2D(double, double);
	Linear_2D(double, double, int, int);
	~Linear_2D();

	const int				n_x, n_y, nx, n_omega;
	const double			x_0, y_0;
	const complex<double>	I;
	double					omega;
	complex<double>			beta;
	complex<double>			f;
	complex<double>			alpha;
	vector<vector<double>>	u;
	vector<vector<double>>	v;
	vector<double>			dp;
	vector<double>			p_curve;
	vector<double>			u_curve;
	vector<double>			omega_vector;

	void					calculate_dp(double);
	void					calculate_u(double);
	void					calculate_v(double);
	void					calculate_u_v(double);
	void					write_in_file(int, string, vector<double>);
	void					write_in_file(int, int, string, vector<vector<double>>);
	void					resonance_curve_p();
	void					resonance_curve_u();
	double					find_abs_max(vector<double>);
	double					find_abs_max(vector<vector<double>>);
};
