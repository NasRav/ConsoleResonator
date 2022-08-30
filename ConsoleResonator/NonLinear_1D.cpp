#include "NonLinear_1D.h"

NonLinear_1D::NonLinear_1D() :
	Resonator(), n_x(200), dx((L + l) / (n_x - 1)), n_t(2 * n_x), alph(1), dt(dx / c0), n_omega(1000)
{
	u.resize(n_x);
	p.resize(n_x);
	du_dx.resize(n_x);
	x.resize(3);
	for (int i = 0; i < 3; i++)
		x[i].resize(n_x);
	for (int i = 0; i < n_x; i++)
	{
		x[0][i] = i * dx;
		x[1][i] = x[0][i];
		u[i] = 0;
		p[i] = 0;
		du_dx[i] = 0;
	}
}

NonLinear_1D::NonLinear_1D(double Length) :
	Resonator(Length), n_x(200), dx((L + l) / (n_x - 1)), n_t(2 * n_x), alph(1), dt(dx / c0), n_omega(1000)
{
	u.resize(n_x);
	p.resize(n_x);
	du_dx.resize(n_x);
	x.resize(3);
	for (int i = 0; i < 3; i++)
		x[i].resize(n_x);
	for (int i = 0; i < n_x; i++)
	{
		x[0][i] = i * dx;
		x[1][i] = x[0][i];
		u[i] = 0;
		p[i] = 0;
		du_dx[i] = 0;
	}
}

NonLinear_1D::NonLinear_1D(double Length, int num_x) :
	Resonator(Length), n_x(num_x), dx(L / (n_x - 1)), n_t(20), alph(1), dt(0.5 * dx / c0), n_omega(1000)
{
	u.resize(n_x);
	p.resize(n_x);
	du_dx.resize(n_x);
	x.resize(3);
	for (int i = 0; i < 3; i++)
		x[i].resize(n_x);
	for (int i = 0; i < n_x; i++)
	{
		x[0][i] = i * dx;
		x[1][i] = x[0][i];
		u[i] = 0;
		p[i] = 0;
		du_dx[i] = 0;
	}
}

NonLinear_1D::~NonLinear_1D()
{}

double	NonLinear_1D::d_dx(vector<double> f, int index)
{
	if (index == 0)
		return -0.5 * (3 * f[0] - 4 * f[1] + f[2]) / dx;
	else if (index == n_x - 1)
		return 0.5 * (3 * f[n_x - 1] - 4 * f[n_x - 2] + f[n_x - 3]) / dx;
	else
		return 0.5 * (f[index + 1] - f[index - 1]) / dx;
}

void	NonLinear_1D::calculate_x_u_p(double omega)
{
	for (int j = 0; j < n_t; j++)
	{
		for (int i = 0; i < n_x; i++)
			du_dx[i] = d_dx(u, i) / d_dx(x[1], i);
		for (int i = 1; i < n_x - 1; i++)
			x[2][i] = 2 * x[1][i] - x[0][i]
			+ dt * dt * c0 * c0 / pow(d_dx(x[1], i), gamma + 1) * (x[1][i + 1] - 2 * x[1][i] + x[1][i - 1]) / (dx * dx)
			+ alph * dt * dt * (4 / 3 * mu * ro0) * d_dx(du_dx, i);
		x[2][0] = 0;
		x[2][n_x - 1] = L + l * cos(omega * j * dt);
		for (int i = 0; i < n_x; i++)
			u[i] = 0.5 * (3 * x[0][i] - 4 * x[1][i] + x[2][i]) / dt;
		for (int i = 0; i < n_x; i++)
		{
			x[0][i] = x[1][i];
			x[1][i] = x[2][i];
		}
	}
	for (int i = 0; i < n_x; i++)
		p[i] = (p0 / pow(d_dx(x[2], i), gamma) - p0) / (ro0 * c0 * omega * l);
}

void	NonLinear_1D::write_in_file(int n, string name, vector<double> array)
{
	ofstream	f_out(name + "_1D_non_linear_L=" + to_string(L) + ".txt");
	for (int i = 0; i < n; i++)
		f_out << array[i] << endl;
	f_out.close();
}

void	NonLinear_1D::write_in_file(int n, string name, vector<double> array1, vector<double> array2)
{
	ofstream	f_out(name + "_1D_non_linear_L=" + to_string(L) + ".txt");
	for (int i = 0; i < n; i++)
		f_out << array1[i] << ' ' << array2[i] << endl;
	f_out.close();
}

void	NonLinear_1D::resonance_curve_u_p()
{
	double	omega_curve;

	omega_vector.resize(2 * n_omega + 1);
	p_curve.resize(2 * n_omega + 1);
	u_curve.resize(2 * n_omega + 1);
	for (int i = -n_omega; i <= n_omega; i++)
	{
		omega_curve = omega_res * (1 + i * (0.1 / n_omega));
		calculate_x_u_p(omega_curve);
		omega_vector[i + n_omega] = omega_curve / omega_res;
		p_curve[i + n_omega] = find_abs_max(p);
		u_curve[i + n_omega] = find_abs_max(u);
	}
}

double	NonLinear_1D::find_abs_max(vector<double> array)
{
	double	max = abs(*array.begin());
	auto	iter = array.begin() + 1;

	while (iter != array.end())
	{
		if (abs(*iter) > max)
			max = abs(*iter);
		++iter;
	}
	return max;
}
