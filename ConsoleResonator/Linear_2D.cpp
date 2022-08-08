#include "Linear_2D.h"

Linear_2D::Linear_2D() :
	n_x(200), n_y(100), I((0, 1))
{
	if (0.5 * L > 700 * delta)
		exit(0);
}

Linear_2D::Linear_2D(int num_x, int num_y) :
	n_x(num_x), n_y(num_y), I((0, 1))
{
	if (0.5 * L > 700 * delta)
		exit(0);
}

Linear_2D::Linear_2D(double Length, double Height) :
	n_x(200), n_y(100), I((0, 1)), Resonator(Length, Height)
{
	if (0.5 * L > 700 * delta)
		exit(0);
}

Linear_2D::Linear_2D(double Length, double Height, int num_x, int num_y) :
	n_x(num_x), n_y(num_y), I((0, 1)), Resonator(Length, Height)
{
	if (0.5 * L > 700 * delta)
		exit(0);
}

Linear_2D::~Linear_2D()
{}
