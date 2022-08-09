#include "Linear_2D.h"

using namespace std;

int main()
{
	Linear_2D	res2d(0.01, 0.01 * 0.03, 6, 6);

	res2d.calculate_dp(res2d.omega_res);
	res2d.calculate_u_v(res2d.omega_res);
	res2d.write_in_file(2 * res2d.nx, "dp", res2d.dp);
	res2d.write_in_file(2 * res2d.nx, res2d.n_y, "u", res2d.u);
	res2d.write_in_file(2 * res2d.nx, res2d.n_y, "v", res2d.v);
	res2d.resonance_curve();
	res2d.write_in_file(2, 2001, "res_cur", res2d.curve);
	return 0;
}
