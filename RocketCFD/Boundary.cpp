#include "Boundary.h"



Boundary::Boundary()
{
}

Boundary::Boundary(BOUNDARY_CONDITION c)
{
	condition = c;
}


Boundary::~Boundary()
{
}

void Boundary::setCondition(BOUNDARY_CONDITION c)
{
	condition = c;
}

void Boundary::update(Face)
{
}
