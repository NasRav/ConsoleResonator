#include "Linear_2D.h"

using namespace std;

int main()
{
	Linear_2D	res2d(1.06, 0.024, 100, 100);

//	res2d.calculate_dp(res2d.omega_res);
//	res2d.calculate_u_v(res2d.omega_res);
//	res2d.write_in_file(2 * res2d.nx + 1, "dp", res2d.dp);
//	res2d.write_in_file(2 * res2d.nx + 1, res2d.n_y, "u", res2d.u);
//	res2d.write_in_file(2 * res2d.nx + 1, res2d.n_y, "v", res2d.v);
	res2d.resonance_curve();
	res2d.write_in_file(3, 2 * res2d.n_omega + 1, "res_cur", res2d.curve);
	return 0;
}
