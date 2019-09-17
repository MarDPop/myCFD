#include "cell.h"
#include <iostream>

Cell::Cell()
{
}

Cell::Cell(const Cell & B)
{
	faces = B.faces;
	neighbors = new Cell * [faces.size()];
	for (size_t i = 0; i < faces.size(); i++) {
		neighbors[i] = B.neighbors[i];
		dx[i] = B.dx[i];
	}
	center = B.center;

}

Cell::~Cell()
{
	delete[] neighbors;
	delete[] dx;
}

/*
void Cell::setFaces(Face ** arr, size_t n)
{
	faces = arr;
	N = n;
	// Faces should be checked for validity (no free vertices, no intersections)
	neighbors = new Cell*[n];
	dx = new Vec3[n];
}
*/

void Cell::collectNeighbors() {
	neighbors = new Cell*[faces.size()];
	dx = new Vec3[faces.size()];
	for (size_t i = 0; i < faces.size(); i++) {
		if (faces[i]->cells[0] == NULL || faces[i]->cells[1] == NULL) {
			neighbors[i] = NULL;
			dx[i] = faces[i]->center - center;
		} else {
			if (faces[i]->cells[0] == this) {
				neighbors[i] = faces[i]->cells[1];
				dx[i] = neighbors[i]->center - center;
			}
			else {
				neighbors[i] = faces[i]->cells[0];
				dx[i] = neighbors[i]->center - center;
			}
		}
	}
}

Face::Face(){}

Face::Face(const Face & u)
{
	N = u.N;
	vertices = new Vec3*[N];
	//std::cout << "Face created \n";
	for (size_t i = 0; i < N; i++) {
		vertices[i] = u[i];
	}
	axial = u.axial;
	center = u.center;
	normal = u.normal;
	area = u.area;
	cells[0] = u.cells[0];
	cells[1] = u.cells[1];
}

Face::Face(Vec3 * v1, Vec3 * v2, bool axi)
{
	vertices = new Vec3*[2];
	vertices[0] = v1;
	vertices[1] = v2;
	N = 2;
	axial = axi;
	center = *v1 + *v2;
	center /= 2;
	double dx = (*v2)[0] - (*v1)[0];
	double dy = (*v2)[1] - (*v1)[1];
	area = sqrt(dx*dx + dy*dy);
	if (axi) {
		area *= 2 * pi * center[1];
		double r1 = (*v1)[1];
		double r2 = (*v2)[1];
		normal[0] = pi*(r1*r1 - r2*r2);	 // "-dy" area
		normal[1] = sqrt(area * area - normal[0] * normal[0]); // Should be "cos" or "dx" of projection but double check
		normal[2] = 0;
	}
	else {
		normal[0] = -dy;
		normal[1] = dx;
		normal[2] = 0;
	}
}

Face::Face(Vec3 * v1, Vec3 * v2, Vec3 * v3)
{
	vertices = new Vec3*[3];
	vertices[0] = v1;
	vertices[1] = v2;
	vertices[2] = v3;
	N = 3;
	center = *v1 + *v2 + *v3;
	center /= 3;
	Vec3 edge1 = *v1 - *v2;
	Vec3 edge2 = *v3 - *v2;
	normal = edge2^edge1;
	normal /= 2;
	area = ~normal;
}

Face::Face(Vec3 ** v, size_t n)
{
	vertices = v;
	N = n;
}

Face::~Face()
{
	std::cout << "Face deleted \n";
	delete[] vertices; // ????
}

Vec3 ** Face::getVertices()
{
	return vertices;
}

size_t Face::getN()
{
	return N;
}

Vec3 * Face::operator[](const size_t & i)
{
	return vertices[i];
}

Vec3 * Face::operator[](const size_t & i) const
{
	return vertices[i];
}

Vec3 Face::operator()(const size_t & i)
{
	return *vertices[i];
}

Vec3 Face::operator*(const Vec3 & v)
{
	return normal*v;
}
