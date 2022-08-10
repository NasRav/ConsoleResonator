#include "Linear_1D.h"

Linear_1D::Linear_1D() :
	Resonator(), n_x(200), n_t(8), n_omega(1000)
{
	u.resize(n_x);
	for (int i = 0; i < n_x; i++)
		u[i].resize(n_t);
	dp.resize(n_x);
	for (int i = 0; i < n_x; i++)
		dp[i].resize(n_t);
}

Linear_1D::Linear_1D(double Length) :
	Resonator(Length, 0), n_x(200), n_t(8), n_omega(1000)
{
	u.resize(n_x);
	for (int i = 0; i < n_x; i++)
		u[i].resize(n_t);
	dp.resize(n_x);
	for (int i = 0; i < n_x; i++)
		dp[i].resize(n_t);
}

Linear_1D::Linear_1D(double Length, int num_x) :
	Resonator(Length, 0), n_x(num_x), n_t(8), n_omega(1000)
{
	u.resize(n_x);
	for (int i = 0; i < n_x; i++)
		u[i].resize(n_t);
	dp.resize(n_x);
	for (int i = 0; i < n_x; i++)
		dp[i].resize(n_t);
}

Linear_1D::~Linear_1D()
{}

void	Linear_1D::calculate_dp(double omega)
{
	k = omega / c0;
	for (int j = 0; j < n_t; j++)
	{
		for (int i = 0; i < n_x; i++)
		{
			dp[i][j] = sin(2 * PI * j / 8 - k * i * L / (n_x - 1)) -
				cos(2 * PI * j / 8 - k * L) * cos(k * i * L / (n_x - 1)) / sin(k * L);
		}
	}
}

void	Linear_1D::calculate_u(double omega)
{
	k = omega / c0;
	for (int j = 0; j < n_t; j++)
	{
		for (int i = 0; i < n_x; i++)
		{
			u[i][j] = sin(2 * PI * j / 8 - k * i * L / (n_x - 1)) -
				sin(2 * PI * j / 8 - k * L) * sin(k * i * L / (n_x - 1)) / sin(k * L);
		}
	}
}

void	Linear_1D::write_in_file(int n, string name, vector<double> array)
{
	ofstream	f_out(name + "_2D_linear_L=" + to_string(L) + "_H=" + to_string(H) + ".txt");
	for (int i = 0; i < n; i++)
		f_out << array[i] << '\n';
	f_out.close();
}

void	Linear_1D::write_in_file(int n, int m, string name, vector<vector<double>> array)
{
	ofstream	f_out(name + "_1D_linear_L=" + to_string(L) + ".txt");
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
			f_out << array[i][j] << '\t';
		f_out << endl;
	}
	f_out.close();
}

void	Linear_1D::resonance_curve_p()
{
	double	omega_curve;

	omega_vector.resize(2 * n_omega + 1);
	p_curve.resize(2 * n_omega + 1);
	for (int i = -n_omega; i <= n_omega; i++)
	{
		omega_curve = omega_res * (1 + i * (0.1 / n_omega));
		calculate_dp(omega_curve);
		omega_vector[i + n_omega] = omega_curve / omega_res;
		p_curve[i + n_omega] = find_abs_max(dp);
	}
}

void	Linear_1D::resonance_curve_u()
{
	double	omega_curve;

	omega_vector.resize(2 * n_omega + 1);
	u_curve.resize(2 * n_omega + 1);
	for (int i = -n_omega; i <= n_omega; i++)
	{
		omega_curve = omega_res * (1 + i * (0.1 / n_omega));
		calculate_u(omega_curve);
		omega_vector[i + n_omega] = omega_curve / omega_res;
		u_curve[i + n_omega] = find_abs_max(u);
	}
}

double	Linear_1D::find_abs_max(vector<vector<double>> array)
{
	double	max = abs(array[0][0]);

	for (int j = 0; j < n_t; j++)
	{
		for (int i = 0; i < n_x; i++)
		{
			if (abs(array[0][0]) > max)
				max = abs(array[i][j]);
		}
	}
	return max;
}
