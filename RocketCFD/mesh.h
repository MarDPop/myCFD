#pragma once
#include <vector>
#include "meshSettings.h"
#include "cell.h"
#include "Boundary.h"
#include "MeshHelper.h"
#include "shapes/Surf2D.h"
#include "util/fmath.h"
#include "Structure.h"
#include <memory>


using namespace fmath;

class Mesh
{
	std::vector<Vec3 *> vertices; // can't be reference
	std::vector<Face *> faces; // can't be reference
	std::vector<Cell *> cells;
	std::vector<std::unique_ptr<Boundary>> boundaries;
	Structure structure;
	MeshSettings settings;
	MeshHelper helper;
public:
	Mesh();
	~Mesh();
	void setSettings(MeshSettings);
	void importRocketShape(std::vector<double> shape, double dx, size_t rDivisions);
	void setICAir(double, double, Vec3, bool);
	void import2DSurface(Surf2D);
	void printCSV(std::string);

};

