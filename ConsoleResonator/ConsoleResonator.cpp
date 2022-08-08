#include <vector>
#include <fstream>
#include "Resonator.h"
#include "Linear_2D.h"

using namespace std;

int main()
{
	const int				n_x(50), n_y(20);
	const double			ro0(1.225), mu(1.82e-5), T_0(288.15), R_gas(8.31), M_mol(0.029), PI(acos(-1)), gamma(1.4),
							p0(ro0 * R_gas * T_0 / M_mol), c0(sqrt(gamma * p0 / ro0));
	const complex<double>	I(0, 1);
	double					x_0, y_0;

	x_0 = 0.5;
	y_0 = 0.0018;

	double					omega = PI * c0 / (2 * x_0);
	double					delta = sqrt(2 * mu / (omega * ro0));

	if (y_0 > 700 * delta)
		exit(0);

	complex<double>			beta = (I + 1.0) / delta;
	complex<double>			f = tanh(beta * y_0) / (beta * y_0);
	complex<double>			alpha = I * omega / c0 / sqrt(1.0 - f);
	vector<vector<double>>	u(2 * n_x, vector<double>(2 * n_y + 1));
	vector<vector<double>>	v(2 * n_x, vector<double>(2 * n_y + 1));
	vector<double>			dp(2 * n_x);

	ofstream	dp_out("dp_2D_linear.txt");

	for (int i = -n_x; i < n_x; i++)
	{
		dp[i + n_x] = -imag(0.5 / sqrt(1.0 - f) *
			(sinh(alpha * static_cast<double>(x_0 * i / n_x)) / cosh(alpha * x_0) -
				cosh(alpha * static_cast<double>(x_0 * i / n_x)) / sinh(alpha * x_0)) *
			exp(I * 2.0 * PI));
		dp_out << dp[i + n_x] << '\t';
	}
	dp_out.close();

	ofstream	u_out("u_2D_linear.txt");
	ofstream	v_out("v_2D_linear.txt");

	for (int j = -n_y; j <= n_y; j++)
	{
		for (int i = -n_x; i < n_x; i++)
		{
			u[i + n_x][j + n_y] = imag(0.5 / (1.0 - f) *
				(cosh(alpha * static_cast<double>(x_0 * i / n_x)) / cosh(alpha * x_0) -
					sinh(alpha * static_cast<double>(x_0 * i / n_x)) / sinh(alpha * x_0)) *
				(1.0 - cosh(beta * static_cast<double>(y_0 * j / n_y)) / cosh(beta * y_0)) *
				exp(I * 2.0 * PI));
			v[i + n_x][j + n_y] = imag(alpha * y_0 * 0.5 * f / (1.0 - f) *
				(sinh(alpha * static_cast<double>(x_0 * i / n_x)) / cosh(alpha * x_0) -
					cosh(alpha * static_cast<double> (x_0 * i / n_x)) / sinh(alpha * x_0)) *
				(sinh(beta * static_cast<double>(y_0 * j / n_y)) / sinh(beta * y_0) -
					static_cast<double>(j / n_y) / y_0) *
				exp(I * 2.0 * PI));
			u_out << u[i + n_x][j + n_y] << '\t';
			v_out << v[i + n_x][j + n_y] << '\t';
		}	//for i n_x
		u_out << endl;
		v_out << endl;
	}	//for j n_y
	u_out.close();
	v_out.close();

	return 0;
}
