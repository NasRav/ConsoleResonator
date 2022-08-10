#pragma once

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>

using namespace std;

class Resonator
{
public:
	Resonator();
	Resonator(double, double);
	~Resonator();

	const double ro0, mu, T_0, R_gas,
		M_mol, PI, gamma, p0, c0;
	double	L, H, omega_res, delta;
};

