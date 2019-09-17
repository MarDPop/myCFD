#include "rocketShape.h"
#include <math.h>
#include "util/fmath.h"
#include <iostream>

rocketShape::rocketShape()
{
}

rocketShape::rocketShape(double _r_throat, double AreaRatio) : r_throat(_r_throat)
{
	r_combustor = r_throat*1.73;
	r_exit = r_throat*sqrt(AreaRatio);
	current_x = 0;
	current_r = r_combustor;
	current_drdx = 0;
	x.push_back(0.0);
	r.push_back(r_combustor);
	drdx.push_back(0.0);
	l_combustor = r_combustor*0.5;
	addHorizontal(l_combustor);
	double rRound = 0.01;
	double angle = -42 * deg2rad;
	addRound(rRound, angle);
	double dR = rRound - rRound*cos(angle);
	toR(r_throat + dR, angle);
	addRound(rRound, -angle);
	addHorizontal(0.001);
	addRound(rRound, 30 * deg2rad);
	sqrtToR(r_exit);
}

rocketShape::rocketShape(int type) 
{
	
}

rocketShape::~rocketShape()
{
}

void rocketShape::update(int t) {
	drdx.push_back(current_drdx);
	r.push_back(current_r);
	x.push_back(current_x);
	type.push_back(t);
}

void rocketShape::addHorizontal(const double length)
{
	current_x += length;
	current_drdx = 0;
	update(0);
}

void rocketShape::toPoint(double x_new, double r_new)
{
	current_drdx = (x_new - current_x) / (r_new - current_r);
	current_r = r_new;
	current_x = x_new;
	update(1);
}

void rocketShape::toR(const double r_new, double angle)
{
	if (r_new > current_r) {
		if (angle < 0) {
			angle = -angle;
		}
	}

	if (r_new < current_r) {
		if (angle > 0) {
			angle = -angle;
		}
	}
	current_drdx = tan(angle);
	current_x += (r_new - current_r) / current_drdx;
	current_r = r_new;
	update(1); 
}

void rocketShape::deltaX(const double dx, double angle)
{
	current_drdx = tan(angle);
	current_r += dx*current_drdx;
	current_x += dx;
	update(1);
}

void rocketShape::addRound(double radius, double angle)
{
	// positive angle up
	double angle1 = atan(current_drdx);
	double angle2 = angle1 + angle;
	double dr = radius*(cos(angle1) - cos(angle2));
	current_drdx = tan(angle2);
	current_x += radius*abs(sin(angle1) - sin(angle2)); 
	if (angle > 0) {
		current_r += dr;
	} else {
		current_r -= dr;
	}
	update(2);
}

void rocketShape::addPoly3rd(double dx_target, double R_target, double drdx_target )
{
	current_r = R_target;
	current_x += dx_target;
	current_drdx = drdx_target;
	update(4);
}

void rocketShape::addPoly5th(double dx_target, double R_target, double drdx_target, double dAdx_target)
{

}

void rocketShape::sqrtToR(double R)
{
	double b = current_x - current_r*0.5 / current_drdx;
	double a = current_r / sqrt(current_x - b);
	current_x = (R / a)*(R / a) + b;
	current_r = R;
	current_drdx = a*0.5 / sqrt(current_x - b);
	update(3);
}

std::vector<double> rocketShape::evenlySpacedR(double dx) {
	bool flag = true;

	for (int i = 0; i < x.size(); i++) {	    
		std::cout << x.at(i) << ", " << r.at(i) << ", " << drdx.at(i) <<  "\n";
	}

	for (int i = 0; i < type.size(); i++) {
		std::cout << type.at(i) << "\n";
	}
	uint16_t shape = 1;
	current_r = r.at(0);
	current_x = 0;
	std::vector<double> output;
	output.push_back(current_r);
	while (shape < x.size()) {
		current_x += dx;
		if (current_x > x.at(shape)) {
			shape++;
		}
		// double test = calcAt(current_x, shape);
		output.push_back(calcAt(current_x,shape));

	}
	return output;
}

double rocketShape::calcAt(double xtest, int loc)
{
	int shapeFunc = -1;
	if (loc < x.size()) {
		shapeFunc = type.at(loc - 1);
	}
	
	if (shapeFunc == 0) {
		return r.at(loc);
	} else if (shapeFunc == 1) {
		double r0 = r.at(loc - 1);
		double dx = xtest - x.at(loc - 1);
		return r0 + dx*drdx.at(loc);
	} else if  (shapeFunc == 2) {
		double angle1 = atan(drdx.at(loc-1));
		double angle2 = atan(drdx.at(loc));
		double dx = abs(x.at(loc) - x.at(loc - 1));
		double radius = dx / abs( sin(angle2) - sin(angle1));
		double y_center = r.at(loc - 1) - radius*cos(angle1);
		double x_center = x.at(loc - 1) - radius*sin(angle1);
		dx = xtest - x_center;
		//double dy = radius*abs(cos(angle2)-cos(angle1));
		double dy = radius*cos(angle1) - sqrt(radius*radius - dx*dx);
		if (angle2 > angle1) {
			return r.at(loc - 1) + dy;
		} else {
			return r.at(loc - 1) - dy;
		}
	} else if (shapeFunc == 3) {
		// Used because dA/dx is smooth
		double b = x.at(loc - 1) - r.at(loc - 1)*0.5 / drdx.at(loc - 1);
		double a = r.at(loc - 1) / sqrt(x.at(loc - 1) - b);
		return a*sqrt(xtest - b);
	}
	else if (shapeFunc == 4) {
		double dr = r.at(loc) - r.at(loc - 1);
		double dx = x.at(loc) - x.at(loc - 1); // x distance
		// f(x) = ax^3 + bx^2 + cx + d
		// f'(x) = 3ax^2 + 2bx + c
		// f(0) = r @ first point
		// f'(0) = drdx @ first point 
		double b_1 = r.at(loc) - (drdx.at(loc - 1)*dx + r.at(loc - 1));	// sol var 1
		double b_2 = drdx.at(loc) - drdx.at(loc - 1);	// sol var 1
		double a = dx*dx*dx;
		double b = dx*dx;
		double c = 3 * b;
		double d = 2 * dx;
		double det = a*d - b*c;
		double a1 = (d*b_1 - c*b_2)/det;
		double a2 = (-b*b_1 + a*b_2) / det;
		dx = xtest - x.at(loc - 1);
		return r.at(loc - 1) + dx*(drdx.at(loc - 1) + dx*(a1 + dx*a2));
	}
	else if (shapeFunc == 5) {
		 // smoother than 3rd order by smoothing out the dA/dx terms
		double dr = r.at(loc) - r.at(loc - 1);
		double dx = x.at(loc) - x.at(loc - 1); // x distance
		double dAdx1 = 2 * pi*r.at(loc - 1)*drdx.at(loc - 1);
		double dAdx2 = 2 * pi*r.at(loc)*drdx.at(loc);
		int shapeFunc1 = 0;
		if (loc > 1) {
			shapeFunc1 = type.at(loc - 2);
		}
		int shapeFunc2 = 0;
		if (loc < type.size()) {
			shapeFunc2 = type.at(loc);
		}
		double curv1, curv2;
		if (shapeFunc1 < 2) {
			curv1 = 0;
		}
		else if(shapeFunc1 == 2) {
			curv1 = (drdx.at(loc - 1) - drdx.at(loc - 2)) / (x.at(loc - 1) - x.at(loc - 2)); // probably not correct
		}
		else if (shapeFunc1 == 3) {
			curv1 = -drdx.at(loc - 1)*drdx.at(loc - 1) / r.at(loc - 1);
		} if (shapeFunc1 == 4) {

		}
		/*
		Matrix yy = Matrix(4, 4);
		double dx2 = dx*dx;
		A.vals[0][1] = dx2*dx2;
		A.vals[0][0] = dx*A.vals[0][1];
		A.vals[0][2] = dx2*dx;
		A.vals[0][3] = dx2;
		A.vals[1][0] = 5 * A.vals[0][1];
		A.vals[1][1] = 4 * A.vals[0][2];
		A.vals[1][2] = 3 * A.vals[0][3];
		A.vals[1][3] = 2 * dx;
		A.vals[2][0] = 20 * A.vals[0][2];
		A.vals[2][1] = 12 * A.vals[0][3];
		A.vals[2][2] = 6 * dx;
		A.vals[2][3] = 2;
		A.vals[3][0] = 0;
		A.vals[3][1] = 0;
		A.vals[3][2] = 0;
		A.vals[2][3] = 2;
				   */
	}
	else {
		return r.back();
	}
}

