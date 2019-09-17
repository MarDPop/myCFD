#include "Flow.h"
#include <cmath>

const double Flow::GAS_R = 8.3144598; // J/mol K
const double Flow::ATM = 101325; // Pa

Flow::Flow(){}

Flow::~Flow(){}

void Flow::setIdeal(double mw, double gam)
{
	_mw = mw;
	_gam = gam;
	_gam1 = (gam + 1) / 2;
	_gam2 = (gam - 1) / 2;
	_R = Flow::GAS_R / mw;
	_cp = gam*_R / (gam - 1);
	_cv = _cp / gam;
}



double Flow::chokePressureRatio(double gamma) {
	return pow(2 / (gamma + 1), gamma / (gamma - 1));
}

double Flow::shockStaticPressureRatio(double mach, double gamma) {
	return (2 * gamma*mach*mach - gamma + 1) / (gamma + 1);
}

double Flow::shockDensityRatio(double mach, double gamma) {
	mach *= mach;
	return (gamma + 1)*mach / ((gamma - 1)*mach + 2);
}

double Flow::shockTotalPressureRatio(double mach, double gamma) {
	return pow(pow(shockDensityRatio(mach, gamma), gamma) / shockStaticPressureRatio(mach, gamma), 1 / (gamma - 1));
}

double Flow::shockMachRatio(double mach, double gamma) {
	double g = gamma - 1;
	mach *= mach;
	return sqrt((g*mach+2) / (2*gamma*mach-g));
}

double Flow::shockTemperatureRatio(double mach, double gamma) {
	return shockStaticPressureRatio(mach, gamma) / shockDensityRatio(mach, gamma);
}

double * Flow::shockRatios(double M2, double gamma) {
	// https://www.grc.nasa.gov/WWW/K-12/airplane/normal.html
	double * out = new double[3];
	out[1] = gamma - 1;
	out[2] = gamma + 1;
	out[0] = (2 * gamma*M2 - out[1]) / out[2];
	out[1] = out[2] * M2 / (out[1] *M2 + 2);
	out[2] = out[0] / out[1]; // M1 ^2
	return out;
}

double Flow::shockPRatioHelp(double p1p0, double rho1rho0, double gamma) {
	// to be used with shockRatios() to get total pressure ratio
	return pow(pow(rho1rho0, gamma) / p1p0, 1 / (gamma - 1));
}

double Flow::isentropicAreaRatio(double mach, double gamma) {
	double g1 = (gamma + 1)/2;
	double g2 = gamma - 1;
	double g3 = g1 / g2;
	return pow(1+g2*0.5*mach*mach,g3)/(mach*pow(g1, g3));
}

double Flow::beta(double mach, double gamma) {
	return 1 / (1+(gamma-1)*mach*mach*0.5);
}

double * Flow::isentropicRatios(double mach, double gamma) {
	// https://www.grc.nasa.gov/WWW/K-12/airplane/isentrop.html
	// array [ temperature, pressure, Area ]
	double * out = new double[3];
	out[2] = gamma - 1;
	out[0] = 1+out[2]*0.5*mach*mach; // not inversed!
	out[1] = pow(out[0], -gamma / out[2]);
	double g = (gamma + 1) / 2;
	out[2] = pow(out[0] / g, g / out[2]) / mach;
	return out;
}

double Flow::chockedMassFlux(double totalPressure, double totalTemperature, double gamma, double mw) {
	double g = gamma + 1;
	return totalPressure*sqrt(gamma*pow(2 / g, g / (gamma - 1))/(Flow::GAS_R/mw*totalTemperature));
}

double Flow::chockedMassFlux(double totalPressure, double gamma, double soundspeed) {
	double g = (gamma + 1) / 2;
	return totalPressure*gamma*pow(g, g / (1-gamma)) / soundspeed;
}

double Flow::pressureFluxRatio(double mach, double gamma) {
	mach *= mach;
	double g1 = (gamma + 1) / 2;
	double g2 = (gamma - 1) / 2;
	g1 = pow(g1, g1 / g2)*mach;
	return 1 / sqrt(g1 + g2*mach*g1);
}

double Flow::exitMachInNozzleShock(double gamma, double throatArea, double exitArea, double p_e, double p_total_throat) {
	// remember ptotal1*Astar1 = ptotal2*Astar2
	// then Aexit/Athroat * pexit/ptotalthroat = Aexit/Astarexit * pexit/ptotalexit
	// note (gamma+1)/(2*(gamma-1)) - gamma/(gamma-1) = - 1/2
	// solve quadratic equation for only positive solutions a = (gamma-1)/2 b = 1 c = (1/areaRatio)^2 / lambda  where lambda = ((gamma+1)/2)^[(gamma+1)/(2*(gamma-1))]
	double g1 = (gamma + 1) / 2;
	double g2 = gamma - 1;
	double c = throatArea* p_total_throat / (exitArea* p_e*pow(g1, g1/g2));
	c *= c;
	return sqrt((-1+sqrt(1+2*g2*c))/g2);
}

double Flow::areaInNozzleShock(double gamma, double throatArea, double exitArea, double p_total_throat, double p_e)
{
	double g1 = (gamma + 1) / 2;
	double g2 = gamma - 1;
	double c = throatArea* p_total_throat / (exitArea* p_e*pow(g1, g1 / g2));
	c *= c;
	double Me = sqrt((-1 + sqrt(1 + 2*c*g2)) / g2);
	double p0e = p_e*pow(1 + 0.5*g2*Me*Me, gamma / g2);
	// secant method
	p0e /= p_total_throat;
	double Ms = 1 / Me;
	Ms *= Ms;
	double dM = 0.1;
	g1 *= 2;
	double f = p0e - pow(pow(g1*Ms / (Ms*g2 + 2), gamma)*g1 / (2*gamma*Ms - g2), 1/g2);
	Ms += dM;
	for (int i = 0; i < 40; i++) {
		double fold = f;
		f = p0e - pow(pow(g1*Ms / (Ms*g2 + 2), gamma)*g1 / (2 * gamma*Ms - g2), 1 / g2);
		dM *= f/(fold-f);
		Ms += dM;
		if (abs(dM) < 1e-5) {
			g1 *= 0.5;
			return pow((1 + g2*0.5*Ms)/g1,g1/g2) / sqrt(Ms);
		}
	}
	return 0.;
}

double Flow::machFromAreaRatio2(double areaRatio, double gamma) {
	double g1 = (gamma + 1) / 2;
	double g3 = g1 / --gamma;
	// other constants
	double ex = g3 - 1;
	gamma *= 0.5;
	// newton method
	//initialize Me guess
	double Me = sqrt(areaRatio);
	// predivide area Ratio to validate function 
	areaRatio *= pow(g1, g3); 
	// if gamma = 1.4 g3 = 3 and ex = 2;
	for (int i = 0; i < 40; i++) {
		double beta = 1 + gamma*Me*Me;
		double f = pow(beta, g3) / Me;
		double dM = (areaRatio - f) / (g1*pow(beta, ex) - f / Me);
		Me += dM;
		if (abs(dM) < 1e-5) {
			return Me;
		}
	}
	return 0.;
}

double Flow::prandtlMeyerAngle(double mach, double gamma)
{
	gamma = sqrt(gamma + 1 / (gamma - 1));
	mach = sqrt(mach*mach - 1);
	return gamma*atan(mach/gamma)-atan(mach);
}

double Flow::machAngle(double mach)
{
	return asin(1/mach);
}

double Flow::vanDerWaalsPressure(double moles, double volume, double temperature, double a, double b)
{
	return moles*(GAS_R*temperature/(volume - moles*b)-a*moles/(volume*volume)) ;
}

double Flow::vanDerWaalsPressure(double molarVolume, double temperature, double a, double b)
{
	return GAS_R*temperature / (molarVolume - b) - a/ (molarVolume*molarVolume);
}

double Flow::vanDerWaalsGasMixA(double * moleFraction, double * a, int n)
{
	double sum = 0;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			sum += moleFraction[i] * moleFraction[j] * sqrt(a[i] * a[j]);
	return sum;
}

double Flow::vanDerWaalsGasMixB(double * moleFraction, double * b, int n)
{
	double sum = 0;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			sum += moleFraction[i] * moleFraction[j] * sqrt(b[i] * b[j]);
	return sum;
}

double Flow::sutherlandViscosity(double temperature, double viscosity_ref, double sutherland_C, double temp_ref ) {
	return  viscosity_ref*((0.555*temp_ref + sutherland_C) / (0.555*temperature + sutherland_C)) * pow(temperature / temp_ref, 1.5);
}

double Flow::machFromAreaRatio(double areaRatio, double gamma) {
	// Generally Faster
	double g1 = (gamma + 1) *0.5;
	double ex = g1 / --gamma;
	// other constants
	gamma *= 0.5;
	// newton method
	//initialize Me guess
	double Me = sqrt(areaRatio);
	// predivide area Ratio to validate function 
	areaRatio *= pow(g1, ex);
	double f = areaRatio - pow(1 + gamma*Me*Me, ex) / Me;
	double dM = 0.01;
	Me += dM;
	for (int i = 0; i < 40; i++) {
		double fold = f;
		f = areaRatio - pow(1 + gamma*Me*Me, ex) / Me;
		dM *= f/(fold - f);
		Me += dM;
		if (abs(dM) < 1e-5) {
			return Me;
		}
	}
	return 0.;
}
