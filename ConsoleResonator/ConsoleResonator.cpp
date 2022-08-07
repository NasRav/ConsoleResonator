#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

int main()
{
	const int				n_t(8), n_x(50), n_y(20);
	const double			ro0(1.225), mu(1.82e-5), T_0(288.15), R_gas(8.31), M(0.029), PI(acos(-1)), gamma(1.4),
							p0(ro0 * R_gas * T_0 / M), c0(sqrt(gamma * p0 / ro0));
	const complex<double>	I(0, 1);
	double					x_0, y_0;

//	cout << "x_0 [m] = ";
//	cin >> x_0;
//	cout << "y_0 [m] = ";
//	cin >> y_0;
	x_0 = 0.5;
//	y_0 = 0.0018;

	double					omega = 2 * PI * c0 / (4 * x_0);
	double					delta = sqrt(2 * mu / (omega * ro0));
	y_0 = 200 * delta;
	if (y_0 > 700 * delta)
		exit(0);
	complex<double>			beta = (I + 1.0) / delta;
	complex<double>			f = tanh(beta * y_0) / (beta * y_0);
	complex<double>			alpha = I * omega / c0 / sqrt(1.0 - f);
	complex<double>			u_complex;
	complex<double>			v_complex;
	complex<double>			dp_complex;
	vector<vector<double>>	u(2 * n_x, vector<double>(2 * n_y));
	vector<vector<double>>	v(2 * n_x, vector<double>(2 * n_y));
	vector<double>			dp(2 * n_x);

	for (int k = 0; k < n_t; k++)
	{
		ofstream	u_out("u_2D_linear_" + to_string(k) + ".txt");
		ofstream	v_out("v_2D_linear_" + to_string(k) + ".txt");
		ofstream	dp_out("dp_2D_linear_" + to_string(k) + ".txt");

		for (int i = -n_x; i < n_x; i++)
		{
			dp_complex = 0.5 / sqrt(1.0 - f) *
				(sinh(alpha * static_cast<complex<double>>(x_0 * i / n_x)) / cosh(alpha * x_0) -
					cosh(alpha * static_cast<complex<double>>(x_0 * i / n_x)) / sinh(alpha * x_0)) *
				exp(I * 2.0 * PI * static_cast<complex<double>>(k / n_t));
			dp[i + n_x] = -imag(dp_complex);
			dp_out << dp[i + n_x] << ' ';
		}
		dp_out.close();

		for (int j = -n_y; j < n_y; j++)
		{
			for (int i = -n_x; i < n_x; i++)
			{
				u_complex = 0.5 / (1.0 - f) *
					(cosh(alpha * static_cast<complex<double>>(x_0 * i / n_x)) / cosh(alpha * x_0) -
						sinh(alpha * static_cast<complex<double>>(x_0 * i / n_x)) / sinh(alpha * x_0)) *
					(1.0 - cosh(beta * static_cast<complex<double>>(y_0 * j / n_y)) / cosh(beta * y_0)) *
					exp(I * 2.0 * PI * static_cast<complex<double>>(k / n_t));
				v_complex = alpha * y_0 * 0.5 * f / (1.0 - f) *
					(sinh(alpha * static_cast<complex<double>>(x_0 * i / n_x)) / cosh(alpha * x_0) -
						cosh(alpha * static_cast<complex<double>> (x_0 * i / n_x)) / sinh(alpha * x_0)) *
					(sinh(beta * static_cast<complex<double>>(y_0 * j / n_y)) / sinh(beta * y_0) -
						static_cast<complex<double>>(j / n_y) / y_0) *
					exp(I * 2.0 * PI * static_cast<complex<double>>(k / n_t));
				u[i + n_x][j + n_y] = imag(u_complex);
				v[i + n_x][j + n_y] = imag(v_complex);
				u_out << u[i + n_x][j + n_y] << ' ';
				v_out << v[i + n_x][j + n_y] << ' ';
			}	//for i n_x
			u_out << endl;
			v_out << endl;
		}	//for j n_y
		u_out.close();
		v_out.close();
	}	//for k n_t

	return 0;
}
