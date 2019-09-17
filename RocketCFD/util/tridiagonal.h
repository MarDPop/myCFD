#pragma once
#include "fmath.h"

using namespace fmath;

struct Tridiagonal
{
	double A,B,C;
	bool preserve, constant;
	size_t n;

	double* a;
	double* b;
	double* c;

	double* c_;
	double* d_;

	Tridiagonal();
	Tridiagonal(double A, double B, double C);
	~Tridiagonal();

	void solve(const Vector & X, Vector & buffer);

	Vector operator/(const Vector &);
	Vector operator*(const Vector &);
	Matrix toMatrix();
};

