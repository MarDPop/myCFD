#pragma once
#include <vector>
#include "../cell.h"

class Shape
{
public:
	std::vector<Vec3 *> vertices;
	std::vector<Face *> faces;

	Shape();
	~Shape();
};

