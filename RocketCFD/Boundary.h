#pragma once
#include "cell.h"
#include <vector>
#include "shapes/Shape.h"

enum BOUNDARY_CONDITION {
	 FREESTREAM,
	 WALL,
	 MASS_FLOW,
	 SYMMETRY,
	 AXIAL,
	 PRESSURE_OUTLET,
	 INTERFACE,
	 SPECIAL // for polymorphed
};

class Boundary
{
	BOUNDARY_CONDITION condition;
	double * values;
	bool * options;
	Boundary * bc; // for interface
public:
	Shape shape;

	Boundary();
	Boundary(BOUNDARY_CONDITION);
	~Boundary();

	void setCondition(BOUNDARY_CONDITION);
	static void setInterface(Boundary, Boundary);
	void update(Face);
};

