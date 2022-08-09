#include "Linear_2D.h"

Linear_2D::Linear_2D() :
	n_x(200), n_y(100), I(0, 1), x_0(0.5 * L), y_0(0.5 * H), nx(0.5 * n_x)
{
	u.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		u[i].resize(n_y);
	v.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		v[i].resize(n_y);
	dp.resize(2 * nx + 1);
}

Linear_2D::Linear_2D(int num_x, int num_y) :
	n_x(num_x), n_y(num_y), I(0, 1), x_0(0.5 * L), y_0(0.5 * H), nx(0.5 * n_x)
{
	u.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		u[i].resize(n_y);
	v.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		v[i].resize(n_y);
	dp.resize(2 * nx + 1);
}

Linear_2D::Linear_2D(double Length, double Height) :
	n_x(200), n_y(100), I(0, 1), Resonator(Length, Height), x_0(0.5 * L), y_0(0.5 * H), nx(0.5 * n_x)
{
	u.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		u[i].resize(n_y);
	v.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		v[i].resize(n_y);
	dp.resize(2 * nx + 1);
}

Linear_2D::Linear_2D(double Length, double Height, int num_x, int num_y) :
	n_x(num_x), n_y(num_y), I(0, 1), Resonator(Length, Height), x_0(0.5 * L), y_0(0.5 * H), nx(0.5 * n_x)
{
	u.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		u[i].resize(n_y);
	v.resize(2 * nx + 1);
	for (int i = 0; i <= 2 * nx; i++)
		v[i].resize(n_y);
	dp.resize(2 * nx + 1);
}

Linear_2D::~Linear_2D()
{}

void	Linear_2D::calculate_dp(double omega)
{
	delta = sqrt(2 * mu / (omega * ro0));
	if (y_0 > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * y_0) / (beta * y_0);
	alpha = I * omega / c0 / sqrt(1.0 - f);

	for (int i = -nx; i <= nx; i++)
	{
		dp[i + nx] = -imag(0.5 / sqrt(1.0 - f) *
			(sinh(alpha * static_cast<double>(x_0 * i / nx)) / cosh(alpha * x_0) -
				cosh(alpha * static_cast<double>(x_0 * i / nx)) / sinh(alpha * x_0)) *
			exp(I * 2.0 * PI));
	}
}

void	Linear_2D::calculate_u(double omega)
{
	delta = sqrt(2 * mu / (omega * ro0));
	if (y_0 > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * y_0) / (beta * y_0);
	alpha = I * omega / c0 / sqrt(1.0 - f);

	for (int j = 0; j < n_y; j++)
	{
		for (int i = -nx; i <= nx; i++)
		{
			u[i + nx][j] = imag(0.5 / (1.0 - f) *
				(cosh(alpha * static_cast<double>(x_0 * i / nx)) / cosh(alpha * x_0) -
					sinh(alpha * static_cast<double>(x_0 * i / nx)) / sinh(alpha * x_0)) *
				(1.0 - cosh(beta * static_cast<double>(y_0 * j / (n_y - 1))) / cosh(beta * y_0)) *
				exp(I * 2.0 * PI));
		}	//for i nx
	}	//for j ny
}

void	Linear_2D::calculate_v(double omega)
{
	delta = sqrt(2 * mu / (omega * ro0));
	if (y_0 > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * y_0) / (beta * y_0);
	alpha = I * omega / c0 / sqrt(1.0 - f);

	for (int j = 0; j < n_y; j++)
	{
		for (int i = -nx; i <= nx; i++)
		{
			v[i + nx][j] = imag(alpha * y_0 * 0.5 * f / (1.0 - f) *
				(sinh(alpha * static_cast<double>(x_0 * i / nx)) / cosh(alpha * x_0) -
					cosh(alpha * static_cast<double>(x_0 * i / nx)) / sinh(alpha * x_0)) *
				(sinh(beta * static_cast<double>(y_0 * j / (n_y - 1))) / sinh(beta * y_0) -
					static_cast<double>(j / (n_y - 1))) *
				exp(I * 2.0 * PI));
		}	//for i nx
	}	//for j ny
}

void	Linear_2D::calculate_u_v(double omega)
{
	delta = sqrt(2 * mu / (omega * ro0));
	if (y_0 > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * y_0) / (beta * y_0);
	alpha = I * omega / c0 / sqrt(1.0 - f);

	for (int j = 0; j < n_y; j++)
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
		}	//for i nx
	}	//for j ny
}

void	Linear_2D::write_in_file(int n, string name, vector<double> array)
{
	ofstream	f_out(name + "_2D_linear_L=" + to_string(L) + "_H=" + to_string(H) + ".txt");
	for (int i = 0; i <= n; i++)
	{
		f_out << array[i] << '\t';
	}
	f_out.close();
}

void	Linear_2D::write_in_file(int n, int m, string name, vector<vector<double>> array)
{
	ofstream	f_out(name + "_2D_linear_L=" + to_string(L) + "_H=" + to_string(H) + ".txt");
	for (int j = 0; j < m; j++)
	{
		for (int i = 0; i <= n; i++)
		{
			f_out << array[i][j] << '\t';
		}	//for i nx
		f_out << endl;
	}	//for j ny
	f_out.close();
}

void	Linear_2D::resonance_curve()
{
	double	omega_curve;
	curve.resize(3);
//	cout << "curve.size() = " << curve.size() << endl;
	for (int i = 0; i < 3; i++)
	{
		curve[i].resize(2001);
//		cout << "curve.size()[" << i << "] = " << curve[i].size() << endl;
	}
	for (int i = -1000; i <= 1000; i++)
	{
		omega_curve = omega_res * (1 + i * 0.001);
//		cout << "omega_curve["<< i << "] = " << omega_curve << endl;
		calculate_dp(omega_curve);
		calculate_u(omega_curve);
		curve[0][i + 1000] = omega_curve / omega_res;
//		cout << "curve[0][" << i << "] = " << curve[0][i + 1] << endl;
		curve[1][i + 1000] = find_abs_max(dp);
//		cout << "curve[1][" << i << "] = " << curve[1][i + 1] << endl;
		curve[2][i + 1000] = find_abs_max(u);
//		cout << "curve[2][" << i << "] = " << curve[2][i + 1] << endl;
	}
}

double	Linear_2D::find_abs_max(vector<double> array)
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

double	Linear_2D::find_abs_max(vector<vector<double>> array)
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
