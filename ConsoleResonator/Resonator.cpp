#include "Resonator.h"

Resonator::Resonator() :
	L(1), H(0.4), ro0(1.225), mu(1.82e-5), T_0(288.15), R_gas(8.31),
	M_mol(0.029), PI(acos(-1)), gamma(1.4),
	p0(ro0 * R_gas * T_0 / M_mol), c0(sqrt(gamma * p0 / ro0)),
	omega_res(PI * c0 / L), delta(sqrt(2 * mu / (omega_res * ro0)))
{}

Resonator::Resonator(double Length, double Height) :
	L(Length), H(Height), ro0(1.225), mu(1.82e-5), T_0(288.15), R_gas(8.31),
	M_mol(0.029), PI(acos(-1)), gamma(1.4),
	p0(ro0 * R_gas * T_0 / M_mol), c0(sqrt(gamma * p0 / ro0)),
	omega_res(PI * c0 / L), delta(sqrt(2 * mu / (omega_res * ro0)))
{}

Resonator::~Resonator()
{}
