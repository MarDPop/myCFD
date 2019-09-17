#include <stdio.h>
#include <iostream>
#include <fstream>
#include "rocketShape.h"
#include "util/fmath.h"
#include "mesh.h"
#include "util/Flow.h"

int main(int argc, const char* argv[])
{
	bool plot = false;

	/*
	Matrix A = Matrix(3, 3);
	A.vals[0][0] = 6;
	A.vals[0][1] = 2;
	A.vals[0][2] = 3;
	A.vals[1][0] = 4;
	A.vals[1][1] = 5;
	A.vals[1][2] = 6;
	A.vals[2][0] = 7;
	A.vals[2][1] = 8;
	A.vals[2][2] = 9;

	Vector b = Vector(3);
	b.vals[0] = 3;
	b.vals[1] = 4;
	b.vals[2] = 5;

	Vector x = b/A;

	std::cout << x.vals[0] << ","  << x.vals[1] << "," << x.vals[2] << "\n\n";
	*/

	double throat_r = 0.01;
	double area_ratio = 8;

	rocketShape rocket = rocketShape(throat_r, area_ratio);
	double dx = 0.0005;
	std::vector<double> r = rocket.evenlySpacedR(dx);

	Mesh mesh;
	double radial_divisions = 10;
	mesh.importRocketShape(r, dx, radial_divisions);
	double p_amb = 100000;
	double t_amb = 300;
	mesh.setICAir(p_amb, t_amb,Vec3(false),true);

	if (plot) {

		std::ofstream myfile;
		myfile.open("geometry.csv");

		for (size_t i = 0; i < r.size(); i++) {
			myfile << i*dx << "," << r.at(i) << " \n";
		}

		mesh.printCSV("meshFaces.csv");
	}

	// calculate properties
	double MW = 0.031;
	double chamber_pressure = 6.7e6; //Pa
	double chamber_temperature = 1750; // K
	double gamma = 1.2;
	double mass_flux = Flow::chockedMassFlux(chamber_pressure, chamber_temperature, gamma, MW);
	double approx_exit_mach = Flow::machFromAreaRatio(area_ratio, gamma);
	double approx_exit_temp = chamber_temperature * Flow::beta(approx_exit_mach, gamma);
	double approx_exit_velocity = sqrt(gamma*Flow::GAS_R / MW * approx_exit_temp)*approx_exit_mach;
	double thrust_approx = mass_flux * pi*throat_r*throat_r*approx_exit_velocity;
	

	system("PAUSE");
}