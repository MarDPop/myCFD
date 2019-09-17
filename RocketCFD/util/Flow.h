#pragma once
#include <math.h>
class Flow
{
protected:
	double _mw, _cp, _cv, _R, _gam, _mach, _stagnationTemp, _stagnationPres, _pres, _temp, _reynolds, _dynamicVisc, _heatConductivity ;
	double _gam1, _gam2, _M2, _v,_mdot, _area;
	int nFluids;
public:
	Flow();
	~Flow();

	void setIdeal(double mw,double gam);

	static const double GAS_R;
	static const double ATM;

	static double chokePressureRatio(double gamma);
	static double shockStaticPressureRatio(double mach, double gamma);
	static double shockDensityRatio(double mach, double gamma);
	static double shockTotalPressureRatio(double mach, double gamma);
	static double shockMachRatio(double mach, double gamma);
	static double shockTemperatureRatio(double mach, double gamma);
	static double * shockRatios(double M2, double gamma);
	static double shockPRatioHelp(double p1p0, double rho1rho0, double gamma);
	static double isentropicAreaRatio(double mach, double gamma);
	static double beta(double mach, double gamma);
	static double * isentropicRatios(double mach, double gamma);
	static double chockedMassFlux(double totalPressure, double totalTemperature, double gamma, double mw);
	static double chockedMassFlux(double totalPressure, double gamma, double soundspeed);
	static double pressureFluxRatio(double mach, double gamma);
	static double exitMachInNozzleShock(double gamma, double throatArea, double exitArea, double p_e, double p_total_throat);
	static double areaInNozzleShock(double gamma, double throatArea, double exitArea, double p_total_throat, double p_e);
	static double machFromAreaRatio(double areaRatio, double gamma);
	static double machFromAreaRatio2(double areaRatio, double gamma);
	static double prandtlMeyerAngle(double mach, double gamma);
	static double machAngle(double mach);
	static double vanDerWaalsPressure(double moles, double volume, double temperature, double a, double b);

	double vanDerWaalsPressure(double molarVolume, double temperature, double a, double b);

	static double vanDerWaalsGasMixA(double * moleFraction, double * a, int n);
	static double vanDerWaalsGasMixB(double * moleFraction, double * b, int n);
	static double sutherlandViscosity(double temperature, double viscosity_ref, double sutherland_C, double temp_ref = 298.15);
};

