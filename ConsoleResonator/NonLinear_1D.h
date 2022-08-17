#pragma once

#include "Resonator.h"

class NonLinear_1D : public Resonator
{
public:
	NonLinear_1D();
	NonLinear_1D(double);
	NonLinear_1D(double, int);
	~NonLinear_1D();

	const int				n_x, n_t, alph, n_omega;
	const double			dx, dt;
	vector<vector<double>>	x;
	vector<double>			u;
	vector<double>			p;
	vector<double>			du_dx;
	vector<double>			p_curve;
	vector<double>			u_curve;
	vector<double>			omega_vector;

	double					d_dx(vector<double>, int);
	void					calculate_x_u_p(double);
	void					write_in_file(int, string, vector<double>);
	void					resonance_curve_u_p();
	double					find_abs_max(vector<double>);
};

