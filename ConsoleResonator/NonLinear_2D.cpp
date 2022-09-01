#include "NonLinear_2D.h"

NonLinear_2D::NonLinear_2D() :
	n_x(200), n_y(100), I(0, 1), Resonator(), x_0(0.5 * L), y_0(0.5 * H), nx(0.5 * n_x), n_omega(1000)
{
	u.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		u[i].resize(n_y);
	v.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		v[i].resize(n_y);
	dp.resize(2 * nx + 1);
}

NonLinear_2D::NonLinear_2D(int num_x, int num_y) :
	n_x(num_x), n_y(num_y), I(0, 1), Resonator(), x_0(0.5 * L), y_0(0.5 * H), nx(0.5 * n_x), n_omega(1000)
{
	u.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		u[i].resize(n_y);
	v.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		v[i].resize(n_y);
	dp.resize(2 * nx + 1);
}

NonLinear_2D::NonLinear_2D(double Length, double Height) :
	n_x(200), n_y(100), I(0, 1), Resonator(Length, Height), x_0(0.5 * L), y_0(0.5 * H), nx(0.5 * n_x), n_omega(1000)
{
	u.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		u[i].resize(n_y);
	v.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		v[i].resize(n_y);
	dp.resize(2 * nx + 1);
}

NonLinear_2D::NonLinear_2D(double Length, double Height, int num_x, int num_y) :
	n_x(num_x), n_y(num_y), I(0, 1), Resonator(Length, Height), x_0(0.5 * L), y_0(0.5 * H), nx(0.5 * n_x), n_omega(1000)
{
	u.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		u[i].resize(n_y);
	u_0.resize(2 * nx + 1);
	du_0_dx.resize(2 * nx + 1);
	Y_x.resize(n_y);
	Y_y.resize(n_y);
	G.resize(2 * nx + 1);
	v.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		v[i].resize(n_y);
	dp.resize(2 * nx + 1);
}

NonLinear_2D::~NonLinear_2D()
{}

void	NonLinear_2D::calculate_dp(double omega)
{
	delta = sqrt(2 * mu / (omega * ro0));
	if (y_0 > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * y_0) / (beta * y_0);
	alpha = I * omega / c0 / sqrt(1.0 - f);

/*	for (int i = -nx; i <= nx; i++)
	{
		dp[i + nx] = -imag(0.5 / sqrt(1.0 - f) *
			(sinh(alpha * static_cast<double>(x_0 * i / nx)) / cosh(alpha * x_0) -
				cosh(alpha * static_cast<double>(x_0 * i / nx)) / sinh(alpha * x_0)) *
			exp(I * 2.0 * PI));
	}*/
}

void	NonLinear_2D::calculate_u(double omega)
{
	delta = sqrt(2 * mu / (omega * ro0));
	if (y_0 > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * y_0) / (beta * y_0);
	alpha = I * omega / c0 / sqrt(1.0 - f);

	for (int i = -nx; i <= nx; i++)
	{
		u_0[i + nx] = 0.5 / (1.0 - f)
			* (cosh(alpha * static_cast<double>(x_0 * i / nx)) / cosh(alpha * x_0) -
				sinh(alpha * static_cast<double>(x_0 * i / nx)) / sinh(alpha * x_0));
		du_0_dx[i + nx] = conj( 0.5 * alpha / (1.0 - f)
			* (sinh(alpha * static_cast<double>(x_0 * i / nx)) / cosh(alpha * x_0) -
				cosh(alpha * static_cast<double>(x_0 * i / nx)) / sinh(alpha * x_0)) );
		G[i + nx] = x_0 / pow(abs(l * omega), 2) * conj(u_0[i + nx]) * du_0_dx[i + nx];
	}
	for (int j = 0; j < n_y; j++)
	{
		Y_x[j] = (1.0 - cosh(beta * static_cast<double>(y_0 * j / (n_y - 1))) / cosh(beta * y_0));
		Y_y[j] = (static_cast<double>(y_0 * j / (n_y - 1)) / y_0
			- sinh(beta * static_cast<double>(y_0 * j / (n_y - 1))) / sinh(beta * y_0));
	}
	V_0 = 2 * pow(abs(l * omega), 2) / (x_0 * omega);
/*
	for (int j = 0; j < n_y; j++)
	{
		for (int i = -nx; i <= nx; i++)
		{
			u[i + nx][j] = V(G, static_cast<double>(y_0 * j / (n_y - 1))) + V_0 * Re * (0.25 * I * (1.0 - f) * conj(G[i + nx]) * conj(Y_x[j])
				+ 1 / y_0 * (3 * A_3(G) * pow(static_cast<double>(y_0 * j / (n_y - 1)), 2) / pow(y_0, 2) + A_1(G[i + nx]));
		}
	}*/
/*	for (int j = 0; j < n_y; j++)
	{
		for (int i = -nx; i <= nx; i++)
		{
			u[i + nx][j] = imag(0.5 / (1.0 - f) *
				(cosh(alpha * static_cast<double>(x_0 * i / nx)) / cosh(alpha * x_0) -
					sinh(alpha * static_cast<double>(x_0 * i / nx)) / sinh(alpha * x_0)) *
				(1.0 - cosh(beta * static_cast<double>(y_0 * j / (n_y - 1))) / cosh(beta * y_0)) *
				exp(I * 2.0 * PI));
		}
	}*/
}

void	NonLinear_2D::calculate_v(double omega)
{
	delta = sqrt(2 * mu / (omega * ro0));
	if (y_0 > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * y_0) / (beta * y_0);
	alpha = I * omega / c0 / sqrt(1.0 - f);

/*	for (int j = 0; j < n_y; j++)
	{
		for (int i = -nx; i <= nx; i++)
		{
			v[i + nx][j] = imag(alpha * y_0 * 0.5 * f / (1.0 - f) *
				(sinh(alpha * static_cast<double>(x_0 * i / nx)) / cosh(alpha * x_0) -
					cosh(alpha * static_cast<double>(x_0 * i / nx)) / sinh(alpha * x_0)) *
				(sinh(beta * static_cast<double>(y_0 * j / (n_y - 1))) / sinh(beta * y_0) -
					static_cast<double>(j / (n_y - 1))) *
				exp(I * 2.0 * PI));
		}
	}*/
}

void	NonLinear_2D::calculate_u_v(double omega)
{
	delta = sqrt(2 * mu / (omega * ro0));
	if (y_0 > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * y_0) / (beta * y_0);
	alpha = I * omega / c0 / sqrt(1.0 - f);

/*	for (int j = 0; j < n_y; j++)
	{
		for (int i = -nx; i <= nx; i++)
		{
			u[i + nx][j] = imag(0.5 / (1.0 - f) *
				(cosh(alpha * static_cast<double>(x_0 * i / nx)) / cosh(alpha * x_0) -
					sinh(alpha * static_cast<double>(x_0 * i / nx)) / sinh(alpha * x_0)) *
				(1.0 - cosh(beta * static_cast<double>(y_0 * j / (n_y - 1))) / cosh(beta * y_0)) *
				exp(I * 2.0 * PI));
			v[i + nx][j] = imag(alpha * y_0 * 0.5 * f / (1.0 - f) *
				(sinh(alpha * static_cast<double>(x_0 * i / nx)) / cosh(alpha * x_0) -
					cosh(alpha * static_cast<double>(x_0 * i / nx)) / sinh(alpha * x_0)) *
				(sinh(beta * static_cast<double>(y_0 * j / (n_y - 1))) / sinh(beta * y_0) -
					static_cast<double>(j / (n_y - 1))) *
				exp(I * 2.0 * PI));
		}
	}*/
}

void	NonLinear_2D::write_in_file(int n, string name, vector<double> array)
{
	ofstream	f_out(name + "_2D_non_linear_L=" + to_string(L) + "_H=" + to_string(H) + ".txt");
	for (int i = 0; i < n; i++)
		f_out << array[i] << '\n';
	f_out.close();
}

void	NonLinear_2D::write_in_file(int n, string name, vector<double> array1, vector<double> array2)
{
	ofstream	f_out(name + "_2D_non_linear_L=" + to_string(L) + ".txt");
	for (int i = 0; i < n; i++)
		f_out << array1[i] << ' ' << array2[i] << endl;
	f_out.close();
}

void	NonLinear_2D::write_in_file(int n, int m, string name, vector<vector<double>> array)
{
	ofstream	f_out(name + "_2D_non_linear_L=" + to_string(L) + "_H=" + to_string(H) + ".txt");
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i < n; i++)
			f_out << array[i][j] << '\t';
		f_out << endl;
	}
	f_out.close();
}

void	NonLinear_2D::resonance_curve_p()
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

void	NonLinear_2D::resonance_curve_u()
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

double	NonLinear_2D::find_abs_max(vector<double> array)
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

double	NonLinear_2D::find_abs_max(vector<vector<double>> array)
{
	double	max = abs(array[0][0]);

	for (int j = 0; j < n_y; j++)
	{
		for (int i = 0; i <= 2 * nx; i++)
		{
			if (abs(array[0][0]) > max)
				max = abs(array[i][j]);
		}
	}
	return max;
}
