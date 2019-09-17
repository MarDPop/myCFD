#pragma once
#include "meshSettings.h"	
#include <vector>
#include <memory>
#include "util/fmath.h"

using namespace fmath;

class Cell;

struct Face {
	Cell * cells[2];

	bool upwind; // true/false for cell 1

	Vec3 center; // center of face
	Vec3 normal; // magnitude is area , normal is in for cell 1

	double area; // normal magnitude (length if 2D)

	// Flow functions
	double density;
	Vec3 v;
	double e, p;

	// Constructors
	Face();
	Face(const Face &);
	Face(Vec3 * v1, Vec3 * v2, bool axi); // 2D Face 
	Face(Vec3 * v1, Vec3 * v2, Vec3 * v3); // 3D triangular
	Face(Vec3**, size_t); //3D polyhedral 
	~Face();

	Vec3 ** getVertices();
	size_t getN();

	Vec3 * operator[](const size_t &);
	Vec3 * operator[](const size_t &) const;
	Vec3 operator()(const size_t &);
	Vec3 operator*(const Vec3 & v);
protected:
	Vec3 ** vertices; // array of vertex pointers
	size_t N; // number vertices
	bool axial;
};

/*
Cells should extend a Fluid or Solid or Plasma etc 
*/
class Cell
{
		 
	Cell ** neighbors; // use pointers since isn't decided until later
	Vec3 * dx; // delta vertex centers
	// size_t N; // number faces/neighbors

public:

	std::vector<Face *> faces; // memory overhead not the biggest deal could use pointer array since REGARDLESS checks to make sure faces are valid will be necessary in future
	
	Vec3 center;

	double volume;

	double density;
	Vec3 rhov;
	double rhoe;

	Vec3 v;
	double e, p, t;

	Cell();
	Cell(const Cell&);
	~Cell();

	// void setFaces(Face **, size_t);

	void collectNeighbors();

};

