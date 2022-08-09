#include "Linear_2D.h"

using namespace std;

int main()
{
	Linear_2D	res2d(0.01, 0.01 * 0.03, 6, 6);

	res2d.calculate_dp();
	res2d.calculate_u_v();
	return 0;
}
