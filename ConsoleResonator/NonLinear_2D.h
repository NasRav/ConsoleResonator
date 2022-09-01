#pragma once

#include "Resonator.h"

class NonLinear_2D : public Resonator
{
public:
	NonLinear_2D();
	NonLinear_2D(int, int);
	NonLinear_2D(double, double);
	NonLinear_2D(double, double, int, int);
	~NonLinear_2D();

	const int				n_x, n_y, nx, n_omega;
	const double			x_0, y_0;
	const complex<double>	I;
	double					omega;
	double					V_0;
	complex<double>			beta;
	complex<double>			f;
	complex<double>			alpha;
	vector<double>			dp;
	vector<double>			p_curve;
	vector<double>			u_curve;
	vector<double>			omega_vector;
	vector<vector<double>>	u;
	vector<complex<double>>	u_0;
	vector<complex<double>>	Y_x;
	vector<complex<double>>	Y_y;
	vector<complex<double>>	G;
	vector<vector<double>>	v;
	vector<complex<complex<double>>>	du_0_dx;

	void					calculate_dp(double);
	void					calculate_u(double);
	void					calculate_v(double);
	void					calculate_u_v(double);
	void					write_in_file(int, string, vector<double>);
	void					write_in_file(int, string, vector<double>, vector<double>);
	void					write_in_file(int, int, string, vector<vector<double>>);
	void					resonance_curve_p();
	void					resonance_curve_u();
	double					find_abs_max(vector<double>);
	double					find_abs_max(vector<vector<double>>);
};
