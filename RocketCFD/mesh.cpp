#include "mesh.h"
#include "MeshHelperDefault.h"
#include <iostream>
#include <fstream>

Mesh::Mesh()
{
	settings = MeshSettings{};
	settings.meshType = DEFAULT;
	settings.prismLayers = false;
	settings.viscous = false;
}

Mesh::~Mesh()
{
}

void Mesh::setSettings(const MeshSettings s)
{
	settings = s;
	switch (s.meshType) {
		case DEFAULT:
			helper = MeshHelperDefault();
			break;
		case TRIMMED:
			break;
	}
}

void Mesh::importRocketShape(std::vector<double> shape, double dx, size_t rDivisions)
{
	double dr;
	Boundary * bc1, * bc2, * bc3, * bc4;
	bc1 = new Boundary();
	bc2 = new Boundary();
	bc3 = new Boundary();
	bc4 = new Boundary();

	dr = shape.front() /rDivisions;

	// Boundaries
	// 1. Mass Boundary (Left)
	bc1->setCondition(MASS_FLOW);
	bc1->shape.vertices.push_back(new Vec3(0, 0, 0));
	for (size_t i = 1; i <= rDivisions; i++) {
		bc1->shape.vertices.push_back(new Vec3(0, i*dr, 0));
		bc1->shape.faces.push_back(new Face(bc1->shape.vertices.at(i - 1), bc1->shape.vertices.at(i), true));
	}

	// 2. Wall Boundary (Top)
	bc2->setCondition(WALL);
	bc2->shape.vertices.push_back(new Vec3(0, shape.at(0), 0));
	for (size_t i = 1; i < shape.size(); i++) {
		bc2->shape.vertices.push_back(new Vec3(i*dx, shape.at(i), 0));
		bc2->shape.faces.push_back(new Face(bc2->shape.vertices.at(i - 1), bc2->shape.vertices.at(i), true));
	}

	// 3. Rocket Outlet (Right)
	double xend = (shape.size() - 1)*dx;
	dr = shape.back()/ rDivisions;
	bc3->setCondition(INTERFACE);
	bc3->shape.vertices.push_back(new Vec3(xend, 0, 0));
	for (size_t i = 1; i <= rDivisions; i++) {
		bc3->shape.vertices.push_back(new Vec3(xend, i*dr, 0));
		bc3->shape.faces.push_back(new Face(bc3->shape.vertices.at(i - 1), bc3->shape.vertices.at(i), true));
	}

	// 4. Axis (Bottom)
	bc4->setCondition(AXIAL);
	bc4->shape.vertices.push_back(new Vec3(0, 0, 0));
	for (size_t i = 1; i < shape.size(); i++) {
		bc4->shape.vertices.push_back(new Vec3(i*dx, 0, 0));
		bc4->shape.faces.push_back(new Face(bc4->shape.vertices.at(i - 1), bc4->shape.vertices.at(i), true));
	}

	// Predefine left face
	Face ** facesTmp = new Face*[rDivisions];
	for (size_t i = 0; i < rDivisions; i++) {
		facesTmp[i] = bc1->shape.faces.at(i);
		facesTmp[i]->cells[0] = NULL;
	}

	// Loop through Interior
	for (size_t i = 1; i < shape.size(); i++) {
		// get radius of shape and divide to get delta 
		dr = shape.at(i)/ rDivisions;
		// get axis symmetrix face at bottom of shape
		Face * bottom = bc4->shape.faces.at(i-1);
		// set first cell to null for all boundary faces
		bottom->cells[0] = NULL;

		for (size_t j = 1; j <= rDivisions; j++) {
			// create new vertex (always top right since left and bottom face already defined 
			Vec3 * v = new Vec3(i*dx, j*dr, 0);	

			// collect faces
			Face * left = facesTmp[j - 1]; // grab face from left
			// bottom face already defined
			// New faces created
			Face * right = new Face((*bottom)[1], v, true);
			Face * top = new Face((*left)[1], v,true);

			// create new Cell
			Cell * c = new Cell();
			// define center
			c->center(i*dx-dx*0.5,j*dr-dr*0.5,0.0);
			// define volume (frustrum boolean subtract : http://mathworld.wolfram.com/ConicalFrustum.html)
			c->volume = (dx / 3)*(right->area + left->area + pi * ((*(*top)[0])[1] * (*(*top)[0])[1] - (*(*bottom)[0])[1] * (*(*bottom)[1])[1]));

			// assign cell to faces
			bottom->cells[1] = c;
			left->cells[1] = c;
			top->cells[0] = c;
			right->cells[0] = c;

			// assign faces to cell
			c->faces.push_back(bottom);
			c->faces.push_back(left);
			c->faces.push_back(right);
			c->faces.push_back(top);

			// add cell to mesh
			cells.push_back(c);

			// push face pointers
			bottom = top;
			facesTmp[j-1] = right;
		}
		// because bottom was set to top at end bottom is top of final cell (alternatively boundary cells may be "material" for heat conduction
		bottom->cells[1] = NULL;
	}
	// set right boundary to null (or interface)!
	for (size_t i = 0; i < rDivisions; i++) {
		facesTmp[i]->cells[1] = NULL;
	}

	boundaries.push_back(std::unique_ptr<Boundary>(bc1));
	boundaries.push_back(std::unique_ptr<Boundary>(bc2));
	boundaries.push_back(std::unique_ptr<Boundary>(bc3));
	boundaries.push_back(std::unique_ptr<Boundary>(bc4));

	for (Cell * c : cells) {
		c->collectNeighbors();
	}
}

void Mesh::setICAir(double a, double b, Vec3 v, bool pres_tmp)
{
	double pressure, temperature, density, energy, enthalpy;
	double R = 267.7;
	double cv = R / 0.4;
	if (pres_tmp) {
		pressure = a;
		density = a / (R*b);
		energy = cv*b;
		temperature = b;
	}
	else {
		pressure = a;
		energy = b;
		density = pressure*energy*0.4;
		temperature = energy / R;
	}

	for (Cell * c : cells) {
		c->e = energy;
		c->p = pressure;
		c->t = temperature;
		c->density = density;
		c->v = v;
		c->rhoe = density*energy;
		c->rhov = v*density;
	}
}

void Mesh::import2DSurface(Surf2D)
{
}

void Mesh::printCSV(std::string fn)
{
	std::ofstream myfile;
	myfile.open(fn);
	for (Cell * c : cells) {
		
		for (Face * f : c->faces) {
			myfile << (*(*f)[0])[0] << "," << (*(*f)[0])[1] << "," << (*(*f)[1])[0] << "," << (*(*f)[1])[1] << "\n";
		}
		
	}
}
