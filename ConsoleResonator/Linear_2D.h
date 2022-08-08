#pragma once

#include <complex>
#include "Resonator.h"

class Linear_2D : public Resonator
{
public:
	Linear_2D();
	Linear_2D(int, int);
	Linear_2D(double, double);
	Linear_2D(double, double, int, int);
	~Linear_2D();

	int							n_x, n_y;
	const std::complex<double>	I;
private:
};

