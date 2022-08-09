#include "Linear_2D.h"

Linear_2D::Linear_2D() :
	n_x(200), n_y(100), I(0, 1)
{
	if (0.5 * H > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * 0.5 * H) / (beta * 0.5 * H);
	alpha = I * omega / c0 / sqrt(1.0 - f);
	u.resize(n_x);
	for (int i = 0; i < n_x; i++)
		u[i].resize(n_y);
	v.resize(n_x);
	for (int i = 0; i < n_x; i++)
		v[i].resize(n_y);
	dp.resize(n_x);
}

Linear_2D::Linear_2D(int num_x, int num_y) :
	n_x(num_x), n_y(num_y), I(0, 1)
{
	if (0.5 * H > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * 0.5 * H) / (beta * 0.5 * H);
	alpha = I * omega / c0 / sqrt(1.0 - f);
	u.resize(n_x);
	for (int i = 0; i < n_x; i++)
		u[i].resize(n_y);
	v.resize(n_x);
	for (int i = 0; i < n_x; i++)
		v[i].resize(n_y);
	dp.resize(n_x);
}

Linear_2D::Linear_2D(double Length, double Height) :
	n_x(200), n_y(100), I(0, 1), Resonator(Length, Height)
{
	if (0.5 * H > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * 0.5 * H) / (beta * 0.5 * H);
	alpha = I * omega / c0 / sqrt(1.0 - f);
	u.resize(n_x);
	for (int i = 0; i < n_x; i++)
		u[i].resize(n_y);
	v.resize(n_x);
	for (int i = 0; i < n_x; i++)
		v[i].resize(n_y);
	dp.resize(n_x);
}

Linear_2D::Linear_2D(double Length, double Height, int num_x, int num_y) :
	n_x(num_x), n_y(num_y), I(0, 1), Resonator(Length, Height)
{
	if (0.5 * H > 700 * delta)
		exit(0);
	beta = (I + 1.0) / delta;
	f = tanh(beta * 0.5 * H) / (beta * 0.5 * H);
	alpha = I * omega / c0 / sqrt(1.0 - f);
	u.resize(2 * n_x);
	for (int i = 0; i < 2 * n_x; i++)
		u[i].resize(2 * n_y + 1);
	v.resize(2 * n_x);
	for (int i = 0; i < 2 * n_x; i++)
		v[i].resize(2 * n_y + 1);
	dp.resize(2 * n_x);
}

Linear_2D::~Linear_2D()
{}

void	Linear_2D::calculate_dp()
{
	ofstream	dp_out("dp_2D_linear.txt");
	for (int i = -n_x; i < n_x; i++)
	{
		dp[i + n_x] = -imag(0.5 / sqrt(1.0 - f) *
			(sinh(alpha * static_cast<double>(0.5 * L * i / n_x)) / cosh(alpha * 0.5 * L) -
				cosh(alpha * static_cast<double>(0.5 * L * i / n_x)) / sinh(alpha * 0.5 * L)) *
			exp(I * 2.0 * PI));
		dp_out << dp[i + n_x] << '\t';
	}
	dp_out.close();
}

void	Linear_2D::calculate_u()
{
	int	nx = n_x / 2;
	int	ny = n_y / 2;

	ofstream	u_out("u_2D_linear.txt");
	for (int j = -ny; j <= ny; j++)
	{
		for (int i = -nx; i < nx; i++)
		{
			u[i + nx][j + ny] = imag(0.5 / (1.0 - f) *
				(cosh(alpha * static_cast<double>(0.5 * L * i / nx)) / cosh(alpha * 0.5 * L) -
					sinh(alpha * static_cast<double>(0.5 * L * i / nx)) / sinh(alpha * 0.5 * L)) *
				(1.0 - cosh(beta * static_cast<double>(0.5 * H * j / ny)) / cosh(beta * 0.5 * H)) *
				exp(I * 2.0 * PI));
			u_out << u[i + nx][j + ny] << '\t';
		}	//for i nx
		u_out << endl;
	}	//for j ny
	u_out.close();
}

void	Linear_2D::calculate_v()
{
	int	nx = n_x / 2;
	int	ny = n_y / 2;

	ofstream	v_out("v_2D_linear.txt");
	for (int j = -ny; j <= ny; j++)
	{
		for (int i = -nx; i < nx; i++)
		{
			v[i + nx][j + ny] = imag(alpha * 0.5 * H * 0.5 * f / (1.0 - f) *
				(sinh(alpha * static_cast<double>(0.5 * L * i / nx)) / cosh(alpha * 0.5 * L) -
					cosh(alpha * static_cast<double> (0.5 * L * i / nx)) / sinh(alpha * 0.5 * L)) *
				(sinh(beta * static_cast<double>(0.5 * H * j / ny)) / sinh(beta * 0.5 * H) -
					static_cast<double>(j / ny) / 0.5 * H) *
				exp(I * 2.0 * PI));
			v_out << u[i + nx][j + ny] << '\t';
		}	//for i nx
		v_out << endl;
	}	//for j ny
	v_out.close();
}

void	Linear_2D::calculate_u_v()
{
	ofstream	u_out("u_2D_linear.txt");
	ofstream	v_out("v_2D_linear.txt");
	for (int j = -n_y; j <= n_y; j++)
	{
		for (int i = -n_x; i < n_x; i++)
		{
			this->u[i + n_x][j + n_y] = imag(0.5 / (1.0 - f) *
				(cosh(alpha * static_cast<double>(0.5 * L * i / n_x)) / cosh(alpha * 0.5 * L) -
					sinh(alpha * static_cast<double>(0.5 * L * i / n_x)) / sinh(alpha * 0.5 * L)) *
				(1.0 - cosh(beta * static_cast<double>(0.5 * H * j / n_y)) / cosh(beta * 0.5 * H)) *
				exp(I * 2.0 * PI));
			this->v[i + n_x][j + n_y] = imag(alpha * 0.5 * H * 0.5 * f / (1.0 - f) *
				(sinh(alpha * static_cast<double>(0.5 * L * i / n_x)) / cosh(alpha * 0.5 * L) -
					cosh(alpha * static_cast<double> (0.5 * L * i / n_x)) / sinh(alpha * 0.5 * L)) *
				(sinh(beta * static_cast<double>(0.5 * H * j / n_y)) / sinh(beta * 0.5 * H) -
					static_cast<double>(j / n_y) / 0.5 * H) *
				exp(I * 2.0 * PI));
			u_out << u[i + n_x][j + n_y] << '\t';
			v_out << u[i + n_x][j + n_y] << '\t';
		}	//for i nx
		u_out << endl;
		v_out << endl;
	}	//for j ny
	u_out.close();
	v_out.close();
}