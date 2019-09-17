#include "tridiagonal.h"
#include <string>

Tridiagonal::Tridiagonal()
{
}

Tridiagonal::Tridiagonal(double A, double B, double C)
{
	this->A = A;
	this->B = B;
	this->C = C;
	constant = true;
}


Tridiagonal::~Tridiagonal()
{
	delete[] a;
	delete[] b;
	delete[] c;
	delete[] c_;
	delete[] d_;
}

void Tridiagonal::solve(const Vector & X, Vector & buffer)
{
	size_t i;
	if (constant) {
		if (n != X.n){
			n = X.n;
			c_ = new double[n-1];
			d_ = new double[n];
			
		}
		c_[0] = C / B;
		d_[0] = X.vals[0] / B;
		for ( i = 1; i < n; i++) {
			double y = B - A*c_[i - 1];
			c_[i] = C / y;
			d_[i] = (X.vals[i] - A*d_[i - 1]) / y;
		}
		buffer.vals[n] = d_[n];
		for ( i = n - 1; i >= 0; i--) {
			buffer.vals[i] = d_[i] - c_[i] * buffer.vals[i+1];
		}
	}
	else {
		if (preserve) {
			c_[0] = c[0] / b[0];
			d_[0] = X.vals[0] / b[0];
			for ( i = 1; i < n; i++) {
				double y = b[i] - a[i]*c_[i - 1];
				c_[i] = c[i] / y;
				d_[i] = (X.vals[i] - a[i]*b[i - 1]) / y;
			}
			buffer.vals[n] = d_[n];
			for ( i = n - 1; i >= 0; i--) {
				buffer.vals[i] = d_[i] - c_[i] * buffer.vals[i + 1];
			}
		}
		else {
			for ( i = 1; i < n; i++) {
				double y = a[i] / b[i - 1];
				b[i] -= y*c[i];
				X.vals[i] -= y*X.vals[i - 1];
			}
			buffer.vals[n] = X.vals[n] / b[n];
			for ( i = n-1; i >= 0; i--) {
				buffer.vals[i] = (X.vals[i] - c[i]*buffer.vals[i+1])/ b[n];
			}
		}
	}
}

Vector Tridiagonal::operator/(const Vector & b)
{
	Vector x = Vector(b.n);
	preserve = true;
	solve(b, x);
	return x;
}

Vector Tridiagonal::operator*(const Vector & b)
{
	Vector x = Vector(b.n);
	if (constant) {
		x.vals[0] = B*b.vals[0] + C*b.vals[1];
		for (size_t i = 1; i < n - 1; i++) {
			x.vals[i] = A*b.vals[i - 1] + B*b.vals[i] + C*b.vals[i+1];
		}
		x.vals[n] = A*b.vals[n-1] + B*b.vals[n];
	}
	return x;
}

Matrix Tridiagonal::toMatrix()
{
	Matrix a = Matrix(n,n);
	if (constant) {
		a.vals[0][0] = B;
		a.vals[0][1] = C;
		for (size_t i = 1; i < n - 1; i++) {
			a.vals[i][i-1] = A;
			a.vals[i][i] = B;
			a.vals[i][i + 1] = C;
		}
		a.vals[n][n-1] = A;
		a.vals[n][n] = B;
	}
	return a;
}
