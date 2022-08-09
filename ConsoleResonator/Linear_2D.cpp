#include "Linear_2D.h"

Linear_2D::Linear_2D() :
	n_x(200), n_y(100), I(0, 1), x_0(0.5 * L), y_0(0.5 * H), nx(0.5 * n_x)
{
	if (y_0 > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * y_0) / (beta * y_0);
	alpha = I * omega / c0 / sqrt(1.0 - f);
	u.resize(2 * nx);
	for (int i = 0; i < 2 * nx; i++)
		u[i].resize(n_y);
	v.resize(2 * nx);
	for (int i = 0; i < 2 * nx; i++)
		v[i].resize(n_y);
	dp.resize(2 * nx);
}

Linear_2D::Linear_2D(int num_x, int num_y) :
	n_x(num_x), n_y(num_y), I(0, 1), x_0(0.5 * L), y_0(0.5 * H), nx(0.5 * n_x)
{
	if (y_0 > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * y_0) / (beta * y_0);
	alpha = I * omega / c0 / sqrt(1.0 - f);
	u.resize(2 * nx);
	for (int i = 0; i < 2 * nx; i++)
		u[i].resize(n_y);
	v.resize(2 * nx);
	for (int i = 0; i < 2 * nx; i++)
		v[i].resize(n_y);
	dp.resize(2 * nx);
}

Linear_2D::Linear_2D(double Length, double Height) :
	n_x(200), n_y(100), I(0, 1), Resonator(Length, Height), x_0(0.5 * L), y_0(0.5 * H), nx(0.5 * n_x)
{
	if (y_0 > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * y_0) / (beta * y_0);
	alpha = I * omega / c0 / sqrt(1.0 - f);
	u.resize(2 * nx);
	for (int i = 0; i < 2 * nx; i++)
		u[i].resize(n_y);
	v.resize(2 * nx);
	for (int i = 0; i < 2 * nx; i++)
		v[i].resize(n_y);
	dp.resize(2 * nx);
}

Linear_2D::Linear_2D(double Length, double Height, int num_x, int num_y) :
	n_x(num_x), n_y(num_y), I(0, 1), Resonator(Length, Height), x_0(0.5 * L), y_0(0.5 * H), nx(0.5 * n_x)
{
	if (y_0 > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * y_0) / (beta * y_0);
	alpha = I * omega / c0 / sqrt(1.0 - f);
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

void	Linear_2D::calculate_dp()
{
	ofstream	dp_out("dp_2D_linear.txt");
	for (int i = -nx; i <= nx; i++)
	{
		dp[i + nx] = -imag(0.5 / sqrt(1.0 - f) *
			(sinh(alpha * static_cast<double>(x_0 * i / nx)) / cosh(alpha * x_0) -
				cosh(alpha * static_cast<double>(x_0 * i / nx)) / sinh(alpha * x_0)) *
			exp(I * 2.0 * PI));
		dp_out << dp[i + nx] << '\t';
	}
	dp_out.close();
}

void	Linear_2D::calculate_u()
{
	ofstream	u_out("u_2D_linear.txt");
	for (int j = 0; j < n_y; j++)
	{
		for (int i = -nx; i <= nx; i++)
		{
			u[i + nx][j] = imag(0.5 / (1.0 - f) *
				(cosh(alpha * static_cast<double>(x_0 * i / nx)) / cosh(alpha * x_0) -
					sinh(alpha * static_cast<double>(x_0 * i / nx)) / sinh(alpha * x_0)) *
				(1.0 - cosh(beta * static_cast<double>(y_0 * j / (n_y - 1))) / cosh(beta * y_0)) *
				exp(I * 2.0 * PI));
			u_out << u[i + nx][j] << '\t';
		}	//for i nx
		u_out << endl;
	}	//for j ny
	u_out.close();
}

void	Linear_2D::calculate_v()
{
	ofstream	v_out("v_2D_linear.txt");
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
			v_out << v[i + nx][j] << '\t';
		}	//for i nx
		v_out << endl;
	}	//for j ny
	v_out.close();
}

void	Linear_2D::calculate_u_v()
{
	ofstream	u_out("u_2D_linear.txt");
	ofstream	v_out("v_2D_linear.txt");
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
			u_out << u[i + nx][j] << '\t';
			v_out << v[i + nx][j] << '\t';
		}	//for i nx
		u_out << endl;
		v_out << endl;
	}	//for j ny
	u_out.close();
	v_out.close();
}
