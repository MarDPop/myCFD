#pragma once
#include <cmath>
#include <vector>

namespace fmath
{
	const double deg2rad = 1.74532925199433e-2;
	const double rad2deg = 57.295779513082323;
	const double twopi = 6.28318530717959;
	const double pi = 3.14159265358979323846;
	const double halfpi = 1.570796326794897;

	float rsqrt(float x);
	double rsqrt(double x);

	double qcos(double x);
	double qsin(double x);
	double qtan(double x);
	double qasin(double x);
	double qexp2(double x);

	double powf(double a, double b);

	float acbrt(float x0);
	float acbrt1(float x0);
	float acbrt2(float x0);
	float acbrt2a(float x0);
	float acbrt2b(float x0);

	struct Vec3 {
		

		Vec3();
		Vec3(const Vec3&);
		Vec3(bool);
		Vec3(double, double, double);
		Vec3(double[3]);
		
		void operator=(const Vec3 &);
		double& operator[](const size_t &);
		const double & operator[](const size_t & i) const;
		Vec3 operator+(const double &);
		void operator+=(const double &);
		Vec3 operator*(const double &);
		void operator*=(const double &);
		Vec3 operator/(const double &);
		void operator/=(const double &);
		Vec3 operator^(const Vec3 &); // cross
		double operator&(const Vec3 &);	// dot
		Vec3 operator*(const Vec3 &); // element muliply
		Vec3 operator/(const Vec3 &); // element divide
		void operator*=(const Vec3 &); // element muliply
		void operator/=(const Vec3 &); // element divide
		Vec3 operator-(const Vec3 & v);
		Vec3 operator+(const Vec3 & v);
		double operator~(); // norm
		void operator-(); // normalize sef
		Vec3 operator!(); // get inverted vector
		Vec3 operator()(); //get normalized vector
		void operator()(double,double,double); //set vals

		void qnormalize();

		Vec3 qgetNorm();

		double qnorm();

	private:
		double vals[3];
	};

	struct Vec3f {
		float vals[3];
		float operator[](const size_t &);
	};

	struct Matrix;

	struct Vector
	{
		double* vals;
		size_t n;
		Vector();
		Vector(const Vector&);
		Vector(size_t);
		Vector(size_t m, double a);
		Vector(double*, size_t);
		~Vector();
		double operator*(const Vector &); // dot product 
		Vector operator%(const Vector &); // cross product
		double norm(double p);
		double operator^(double p);
		double operator~(); // norm
		void qnorm3();
		Vector operator+(const Vector &);
		Vector operator-(const Vector & v);
		Vector operator+(const double & v);
		Vector operator-(const double & v);
		Vector operator*(const double & a);
		Vector operator*=(const Vector & a);
		Vector operator/(const Matrix & A);	  // Matrix solv Ax = b;
	};

	struct Matrix
	{
		double** vals;
		size_t rows, cols;
		Matrix();
		~Matrix();
		Matrix(double ** v, size_t r, size_t c);
		Matrix(double ** v, size_t r1, size_t c1, size_t r2, size_t c2);
		Matrix(size_t, size_t);
		Matrix(const Matrix&);
		void operator=(const Matrix & A);
		Matrix operator*=(const Matrix & A); // element multiply
		Matrix operator*(const Matrix &); // matrix multiply (dot product if two vectors)
		Vector operator*(const Vector &);
		Matrix operator*(const double &);
		void mult(const double & a);
		Matrix operator+(const Matrix &); // add
		Matrix operator-(const Matrix &); // subtract
		Vector operator/(const Vector &); //Ax = b solve
		Matrix operator!(); // Inv
		Matrix minor(size_t r, size_t c);

		Matrix operator~();
		Matrix operator+(); // transpose

		Matrix operator%(const Matrix &); // element multiply
		Matrix mult(const Matrix &); // square tiled multiply
		double det();
		double det(size_t i, bool * js);

	private:
		void minor(size_t r, size_t c, Matrix * B);
		double minorDet(size_t r, size_t c, Matrix * B);
	};

	

};





