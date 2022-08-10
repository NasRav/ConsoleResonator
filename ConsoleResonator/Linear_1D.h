#pragma once

#include "Resonator.h"

class Linear_1D : public Resonator
{
public:
	Linear_1D();
	Linear_1D(double);
	Linear_1D(double, int);
	~Linear_1D();

	const int				n_x, n_t, n_omega;
	double					k;
	double					omega;
	vector<vector<double>>	dp;
	vector<vector<double>>	u;
	vector<double>			p_curve;
	vector<double>			u_curve;
	vector<double>			omega_vector;

	void					calculate_dp(double);
	void					calculate_u(double);
	void					write_in_file(int, string, vector<double>);
	void					write_in_file(int, int, string, vector<vector<double>>);
	void					resonance_curve_p();
	void					resonance_curve_u();
	double					find_abs_max(vector<vector<double>>);
};

