#pragma once

#include <complex>
#include <vector>
#include <cmath>
#include <fstream>
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

	const int				n_x, n_y, nx;
	const double			x_0, y_0;
	const complex<double>	I;
	complex<double>			beta;
	complex<double>			f;
	complex<double>			alpha;
	vector<vector<double>>	u;
	vector<vector<double>>	v;
	vector<double>			dp;

	void					calculate_dp();
	void					calculate_u();
	void					calculate_v();
	void					calculate_u_v();
	void					write_in_file(int, string, vector<double>);
	void					write_in_file(int, int, string, vector<vector<double>>);
	void					resonance_curve();
private:
};
