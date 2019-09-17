#pragma once
#include <vector>
#include "util/fmath.h"

using namespace fmath;

class rocketShape
{
	std::vector<double> x;
	std::vector<int> type;
	std::vector<double> r;
	std::vector<double> drdx;

	// radius of throat, exit and combustor
	double r_throat, r_exit, r_combustor;

	double l_combustor, l_throat;

	double ang_combustor, ang_nozzle;

	double current_x, current_r, current_drdx;
	
public:
	enum NOZZLE_SHAPE { CONE, PARABOLA, THIRD_ORDER, HYPERBOLIC, EXPONENTIAL, SPLINE, METHOD_OF_CHARACTERSTICS };
	enum curveType {horizontal, angle, linear, round, poly2, poly3 };

	rocketShape();
	rocketShape(double, double);
	rocketShape(int type);
	~rocketShape();

	void addHorizontal(double);
	void toPoint(double, double);
	void toR(double,double);
	void update(int);
	void deltaX(double, double);
	std::vector<double> evenlySpacedR(double dx);
	double calcAt(double x, int shapeFunc);
	void addRound(double, double);
	void addPoly3rd(double, double, double);
	void addPoly5th(double, double,double, double);
	void sqrtToR(double);
};

