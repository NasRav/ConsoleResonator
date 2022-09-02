#include "Linear_2D.h"
#include "Linear_1D.h"
#include "NonLinear_1D.h"
#include "NonLinear_2D.h"

int main()
{
//	Linear_2D		res2d(1.06, 0.024, 100, 100);
//	Linear_1D		res1d(1.06, 100);
	NonLinear_1D	lagrange(1.06, 100);

//	res1d.resonance_curve_p();
//	res1d.write_in_file(2 * res1d.n_omega + 1, "omega_cur", res1d.omega_vector);
//	res1d.write_in_file(2 * res1d.n_omega + 1, "dp_res_cur", res1d.p_curve);
//	res1d.write_in_file(2 * res1d.n_omega + 1, "dp_res_cur", res1d.omega_vector, res1d.p_curve);

//	res2d.resonance_curve_p();
//	res2d.write_in_file(2 * res2d.n_omega + 1, "omega_cur", res2d.omega_vector);
//	res2d.write_in_file(2 * res2d.n_omega + 1, "dp_res_cur", res2d.p_curve);
//	res2d.write_in_file(2 * res2d.n_omega + 1, "dp_res_cur", res2d.omega_vector, res2d.p_curve);

	lagrange.resonance_curve_u_p();
//	lagrange.write_in_file(2 * lagrange.n_omega + 1, "omega_cur", lagrange.omega_vector);
//	lagrange.write_in_file(2 * lagrange.n_omega + 1, "dp_res_cur", lagrange.p_curve);
	lagrange.write_in_file(2 * lagrange.n_omega + 1, "dp_res_cur", lagrange.omega_vector, lagrange.p_curve);

	return 0;
}
