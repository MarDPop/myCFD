#include <immintrin.h>
#define intrin_ZERO(a,n){\
size_t x = 0;\
const size_t inc = 32 / sizeof(*(a));\
for (;x < n-inc;x+=inc)\
    _mm256_storeu_ps((float *)((a)+x),_mm256_setzero_ps());\
	switch (n - x) {\
	case 3:\
	(a)[x] = 0; x++; \
	case 2:\
	_mm_storeu_ps((float *)((a)+x), _mm_setzero_ps()); break; \
	case 1:\
	(a)[x] = 0; \
	break; \
	case 0:\
	break; \
	}; \
}

#include "fmath.h"
#include <math.h>

namespace fmath {

float rsqrt(float x) {
	float xhalf = 0.5f * x;
	int i = *(int*)&x;
	i = 0x5f3759df - (i >> 1);
	x = *(float*)&i;  x = x*(1.5f - (xhalf*x*x));
	return x;
}

double rsqrt(double x)
{
	double y = x;
	double x2 = y * 0.5;
	std::int64_t i = *(std::int64_t *) &y;
	// The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
	i = 0x5fe6eb50c7b537a9 - (i >> 1);
	y = *(double *)&i;
	y = y * (1.5 - (x2 * y * y));   // 1st iteration
	//      y  = y * ( 1.5 - ( x2 * y * y ) );   // 2nd iteration, this can be removed
	return y;
}

double * cross(Vec3 u, Vec3 v) {
	double * out = new double[3];

	return out;
}

double dot3(Vec3 u, Vec3 v) {
	return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

double norm(Vec3 u) {
	return sqrt(dot3(u, u));
}

void multiply(double a, double u[3]) {
	u[0] *= a;
	u[1] *= a;
	u[2] *= a;
}

void normalize(double u[3]) {
	multiply(1 / sqrt(dot3(u, u)), u);
}

double qcos(double x) {
	return halfpi - qsin(x);
}

double qsin(double x)
{
	if (x < 0) {
		return -sin(-x);
	}
	if (x > halfpi) {
		if (x > pi) {
			if (x > twopi) {
				return sin(fmod(x,twopi));
			}
			return -sin(x - pi);
		}
		return sin(pi - x);
	}

	double res = x;
	double x2 = -x*x;
	x *= x2;
	res += x*0.166666666666;

	if (x2 < -0.01) {
		x *= x2;
		res += x*8.3333333333e-3;
		if (x2 < -0.25f) {
			x *= x2;
			res += x*1.984126984126984e-04;
			if (x < -0.8f) {
				x *= x2;
				res += x*2.755731922398589e-06;
			}
		}
	}
	return res;
}

double qtan(double x) {
	double x2 = x*x;
	return x*(1 - 0.095238095238095*x2) / (1 - 0.428571428571429*x2 + 0.009523809523810*(x2*x2));
}

double qasin(double x) {
	bool negate = x < 0;
	x = abs(x);
	double ret = -0.0187293;
	ret *= x;
	ret += 0.0742610;
	ret *= x;
	ret -= 0.2121144;
	ret *= x;
	ret += 1.5707288;
	ret = halfpi - sqrt(1.0 - x)*ret;
	if (negate) {
		return -ret;
	}
	return  ret;
}

double qexp2(double x) {
	x = 1 + x*2.44140625e-04;
	x *= x; x *= x; x *= x; x *= x;
	x *= x; x *= x; x *= x; x *= x;
	x *= x; x *= x; x *= x;
	return x*x;
}


double powf(double a, double b) {
	union {
		double d;
		int x[2];
	} u = { a };
	u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
	u.x[0] = 0;
	return u.d;
}

/* ----------------------------- acbrt ------------------------------ */

/* This is a novel and fast routine for the approximate cube root of
an IEEE float (single precision). It is derived from a similar program
for the approximate square root.
The relative error ranges from 0 to +0.00103.
Caution:
Result for -0 is NaN.
Result for 0 is 1.24e-13.
For denorms it is either within tolerance or gives a result < 2.1e-13.
Gives the correct result (inf) for x = inf.
Gives the correct result (NaN) for x = NaN. */

float acbrt(float x0) {
	union { int ix; float x; };

	x = x0;                      // x can be viewed as int.
	ix = ix / 4 + ix / 16;           // Approximate divide by 3.
	ix = ix + ix / 16;
	ix = ix + ix / 256;
	ix = 0x2a5137a0 + ix;        // Initial guess.
	x = 0.33333333f*(2.0f*x + x0 / (x*x));  // Newton step.
	return x;
}

/* ----------------------------- acbrt1 ----------------------------- */

/* This is acbrt with an additional step of the Newton iteration, for
increased accuracy.
The relative error ranges from 0 to +0.00000116. */

float acbrt1(float x0) {
	union { int ix; float x; };

	x = x0;                      // x can be viewed as int.
	ix = ix / 4 + ix / 16;           // Approximate divide by 3.
	ix = ix + ix / 16;
	ix = ix + ix / 256;
	ix = 0x2a5137a0 + ix;        // Initial guess.
	x = 0.33333333f*(2.0f*x + x0 / (x*x));  // Newton step.
	x = 0.33333333f*(2.0f*x + x0 / (x*x));  // Newton step again.
	return x;
}

/* ----------------------------- acbrt2 ----------------------------- */

/* This is a very approximate but very fast version of acbrt. It is just
two integer instructions (shift right and divide), plus instructions
to load the constant.
The constant 0x2a51067f balances the relative error at +-0.0316. */

float acbrt2(float x0) {
	union { int ix; float x; };

	x = x0;                      // x can be viewed as int.
	ix = 0x2a51067f + ix / 3;      // Initial guess.
	return x;
}

/* ---------------------------- acbrt2a ----------------------------- */

/* This is a very approximate but very fast version of acbrt. It is just
eight integer instructions (shift rights and adds), plus instructions
to load the constant.
1/3 is approximated as 1/4 + 1/16 + 1/64 + 1/256 + ... + 1/65536.
The constant 0x2a511cd0 balances the relative error at +-0.0321. */

float acbrt2a(float x0) {
	union { int ix; float x; };

	x = x0;                      // x can be viewed as int.
	ix = ix / 4 + ix / 16;           // Approximate divide by 3.
	ix = ix + ix / 16;
	ix = ix + ix / 256;
	ix = 0x2a511cd0 + ix;        // Initial guess.
	return x;
}

/* ---------------------------- acbrt2b ----------------------------- */

/* This is a very approximate but very fast version of acbrt. It is just
six integer instructions (shift rights and adds), plus instructions
to load the constant.
1/3 is approximated as 1/4 + 1/16 + 1/64 + 1/256.
The constant 0x2a6497f8 balances the relative error at +-0.151. */

float acbrt2b(float x0) {
	union { int ix; float x; };

	x = x0;                      // x can be viewed as int.
	ix = ix / 4 + ix / 16;           // Approximate divide by 3.
	ix = ix + ix / 16;
	ix = 0x2a6497f8 + ix;        // Initial guess.
	return x;
}


Vec3::Vec3(double u[3])
{
	std::memcpy(vals,u,24);
}

Vec3::Vec3(){}

Vec3::Vec3(const Vec3 & u)
{
	std::memcpy(vals, u.vals, 24);
}

Vec3::Vec3(double x, double y, double z)
{
	vals[0] = x;
	vals[1] = y;
	vals[2] = z;
}

Vec3::Vec3(bool b)
{
	memset(vals, b, 24);
}

void Vec3::operator=(const Vec3 & u)
{
	std::memcpy(vals, u.vals, 24);
}

	double& Vec3::operator[](const size_t & i)
	{
		return vals[i];
	}

	const double& Vec3::operator[](const size_t & i) const
	{
		return vals[i];
	}

	Vec3 Vec3::operator+(const double & a)
	{
		Vec3 u;
		u[0] = vals[0] + a;
		u[1] = vals[1] + a;
		u[2] = vals[2] + a;
		return u;
	}

	void Vec3::operator+=(const double & a)
	{
		vals[0] += a;
		vals[1] += a;
		vals[2] += a;
	}

	Vec3 Vec3::operator*(const double & a)
	{
		Vec3 u;
		u[0] = vals[0] * a;
		u[1] = vals[1] * a;
		u[2] = vals[2] * a;
		return u;
	}

	void Vec3::operator*=(const double & a)
	{
		vals[0] *= a;
		vals[1] *= a;
		vals[2] *= a;
	}

	Vec3 Vec3::operator/(const double & a)
	{
		Vec3 u;
		u[0] = vals[0] / a;
		u[1] = vals[1] / a;
		u[2] = vals[2] / a;
		return u;
	}

	void Vec3::operator/=(const double & a)
	{
		vals[0] /= a;
		vals[1] /= a;
		vals[2] /= a;
	}

	Vec3 Vec3::operator^(const Vec3 & v)
	{
		Vec3 out;
		out[0] = vals[1] * v[2] - vals[2] * v[1];
		out[1] = vals[2] * v[0] - vals[0] * v[2];
		out[2] = vals[0] * v[1] - vals[1] * v[0];
		return out;
	}

	double Vec3::operator&(const Vec3 & v)
	{
		return vals[0] * v[0]+ vals[1] * v[1]+ vals[2] * v[2];
	}

	Vec3 Vec3::operator*(const Vec3 & v)
	{
		Vec3 out;
		out[0] = vals[0] * v[0];
		out[1] = vals[1] * v[1];
		out[2] = vals[2] * v[2];
		return out;
	}

	Vec3 Vec3::operator/(const Vec3 & v)
	{
		Vec3 out;
		out[0] = vals[0] / v[0];
		out[1] = vals[1] / v[1];
		out[2] = vals[2] / v[2];
		return out;
	}

	void Vec3::operator*=(const Vec3 & v)
	{
		vals[0] *= v[0] ;
		vals[1] *= v[1] ;
		vals[2] *= v[2] ;
	}

	void Vec3::operator/=(const Vec3 & v)
	{
		vals[0] /= v[0];
		vals[1] /= v[1];
		vals[2] /= v[2];
	}

	Vec3 Vec3::operator-(const Vec3 & v)
	{
		Vec3 out;
		out[0] = vals[0] - v[0];
		out[1] = vals[1] - v[1];
		out[2] = vals[2] - v[2];
		return out;
	}

	Vec3 Vec3::operator+(const Vec3 & v)
	{
		Vec3 out;
		out[0] = vals[0] + v[0];
		out[1] = vals[1] + v[1];
		out[2] = vals[2] + v[2];
		return out;
	}

	double Vec3::operator~()
	{
		return sqrt(vals[0] * vals[0] + vals[1] * vals[1] + vals[2] * vals[2]);
	}

	void Vec3::operator-()
	{
		vals[0] = -vals[0];
		vals[1] = -vals[1];
		vals[2] = -vals[2];
	}

	Vec3 Vec3::operator()()
	{
		Vec3 out;
		double y = 1 / sqrt(vals[0] * vals[0] + vals[1] * vals[1] + vals[2] * vals[2]);
		out[0] = vals[0] * y;
		out[1] = vals[1] * y;
		out[2] = vals[2] * y;
		return out;
	}

	void Vec3::operator()(double x, double y, double z)
	{
		vals[0] = x;
		vals[1] = y;
		vals[2] = z;
	}

	Vec3 Vec3::operator!()
	{
		double y = 1 / sqrt(vals[0] * vals[0] + vals[1] * vals[1] + vals[2] * vals[2]);
		vals[0] /= y;
		vals[1] /= y;
		vals[2] /= y;
		return true;
	}

	void Vec3::qnormalize() {
		double y = rsqrt(vals[0] * vals[0] + vals[1] * vals[1] + vals[2] * vals[2]);
		vals[0] /= y;
		vals[1] /= y;
		vals[2] /= y;
	}

	Vec3 Vec3::qgetNorm() {
		Vec3 out;
		double y = rsqrt(vals[0] * vals[0] + vals[1] * vals[1] + vals[2] * vals[2]);
		out[0] = vals[0] * y;
		out[1] = vals[1] * y;
		out[2] = vals[2] * y;
		return out;
	}

	double Vec3::qnorm() {
		return rsqrt(vals[0] * vals[0] + vals[1] * vals[1] + vals[2] * vals[2]);
	}

	float Vec3f::operator[](const size_t & i)
	{
		return vals[i];
	}

	//-------------------------------------------------------------------------------------------

	Vector::Vector()
	{
	}

	Vector::Vector(const Vector & v)
	{
		vals = new double[v.n];
		n = v.n;
		std::memcpy(vals, v.vals, 8*n);
	}

	Vector::Vector(size_t m)
	{
		vals = new double[m];
		n = m;
	}

	Vector::Vector(size_t m, double a)
	{
		vals = new double[m];
		std::fill_n(vals, m, a);
		n = m;
	}

	Vector::Vector(double * v, size_t m)
	{
		vals = v;
		n = m;
	}

	Vector::~Vector()
	{
		delete[] vals;
	}

	double Vector::operator*(const Vector & v)
	{
		double sum = 0;
		for (size_t i = 0; i < n; i++)
			sum += v.vals[i] * vals[i];
		return sum;
	}

	Vector Vector::operator%(const Vector & v)
	{
		double* x = new double[3];
		x[0] = vals[1] * v.vals[2] - vals[2] * v.vals[1];
		x[1] = vals[2] * v.vals[0] - vals[0] * v.vals[2];
		x[2] = vals[0] * v.vals[1] - vals[1] * v.vals[0];
		return Vector(x, 3);
	}

	double Vector::norm(double p) {
		double sum = 0;
		for (size_t i = 0; i < n; i++)
			sum += pow(vals[i], p);
		return pow(sum, 1 / p);
	}

	double Vector::operator~() {
		double sum = 0;
		for (size_t i = 0; i < n; i++)
			sum += vals[i] * vals[i];
		return sqrt(sum);
	}

	void Vector::qnorm3()
	{
		double y = vals[0] * vals[0] + vals[1] * vals[1] + vals[2] * vals[2];
		double x2 = y * 0.5;
		std::int64_t i = *(std::int64_t *) &y;
		// The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
		i = 0x5fe6eb50c7b537a9 - (i >> 1);
		y = *(double *)&i;
		y = y * (1.5 - (x2 * y * y));   // 1st iteration
										//      y  = y * ( 1.5 - ( x2 * y * y ) );   // 2nd iteration, this can be removed
		vals[0] *= y;
		vals[1] *= y;
		vals[2] *= y;
	}

	Vector Vector::operator+(const Vector & v)
	{
		Vector u = Vector(n);
		for (size_t i = 0; i < n; i++) {
			u.vals[i] = vals[i] + v.vals[i];
		}
		return u;
	}

	Vector Vector::operator-(const Vector & v)
	{
		Vector u = Vector(n);
		for (size_t i = 0; i < n; i++) {
			u.vals[i] = vals[i] - v.vals[i];
		}
		return u;
	}

	Vector Vector::operator+(const double & a)
	{
		Vector u = Vector(n);
		for (size_t i = 0; i < n; i++) {
			u.vals[i] = vals[i] + a;
		}
		return u;
	}

	Vector Vector::operator-(const double & a)
	{
		Vector u = Vector(n);
		for (size_t i = 0; i < n; i++) {
			u.vals[i] = vals[i] - a;
		}
		return u;
	}

	Vector Vector::operator*(const double & a)
	{
		Vector u = Vector(n);
		for (size_t i = 0; i < n; i++) {
			u.vals[i] = vals[i] * a;
		}
		return u;
	}

	Vector Vector::operator*=(const Vector & a)
	{
		Vector u = Vector(n);
		for (size_t i = 0; i < n; i++) {
			u.vals[i] = vals[i] * a.vals[i];
		}
		return Vector();
	}

	Vector Vector::operator/(const Matrix & A)
	{
		// b.n must equal rows
		double* x = new double[n];

		double maxA, *ptr;
		size_t imax, i, j, k;

		size_t * p = new size_t[n];
		double** a = new double*[n];
		for (i = 0; i < n; i++) {
			p[i] = i;
			a[i] = new double[n];
			std::memcpy(a[i], A.vals[i], A.cols * 8);
		}

		for (i = 0; i < n; ++i) {
			maxA = 0.0;
			imax = i;

			for (k = i; k < n; k++)
				if (abs(a[k][i]) > maxA) {
					maxA = abs(a[k][i]);
					imax = k;
				}

			if (maxA < 1e-64) return false;

			if (imax != i) {
				//pivoting P
				j = p[i];
				p[i] = p[imax];
				p[imax] = j;

				//pivoting rows of A
				ptr = a[i];
				a[i] = a[imax];
				a[imax] = ptr;
			}

			for (j = i + 1; j < n; j++) {
				a[j][i] /= a[i][i];
				for (k = i + 1; k < n; k++)
					a[j][k] -= a[j][i] * a[i][k];
			}
		}

		for (i = 0; i < n; i++) {
			x[i] = vals[p[i]];

			for (k = 0; k < i; k++)
				x[i] -= a[i][k] * x[k];
		}

		for (int i = n - 1; i >= 0; i--) {
			for (k = i + 1; k < n; k++)
				x[i] -= a[i][k] * x[k];

			x[i] /= a[i][i];
		}
		delete[] p;
		for (i = 0; i<n; i++)
			delete[] a[i];
		delete[] a;

		return Vector(x, n);
	}


	//-------------------------------------------------------------------------------------------


	Matrix::Matrix()
	{
	}

	Matrix::~Matrix()
	{
		for (size_t i = 0; i<rows; i++)
			delete[] vals[i];
		delete[] vals;
	}

	Matrix::Matrix(double** v, size_t r, size_t c)
	{
		rows = r;
		cols = c;
		vals = v;
	}

	Matrix::Matrix(double** v, size_t r1, size_t c1, size_t r2, size_t c2)
	{
		//https://stackoverflow.com/questions/29375797/copy-2d-array-using-memcpy
		rows = r2 - r1 + 1;
		cols = c2 - c1 + 1;
		vals = new double*[rows];
		for (int i = 0; i < rows; i++) {
			vals[i] = new double[cols];
			std::memcpy(vals[i], v[i + r1] + c1, cols * 8);
		}
	}

	Matrix::Matrix(size_t r, size_t c)
	{
		rows = r;
		cols = c;
		vals = new double*[r];
		for (size_t i = 0; i < r; i++)
			vals[i] = new double[c];
	}

	Matrix::Matrix(const Matrix & A)
	{
		rows = A.rows;
		cols = A.cols;
		vals = new double*[rows];
		for (size_t i = 0; i < rows; i++)
			vals[i] = new double[cols];
		memcpy(vals, A.vals, sizeof(vals));
	}

	void Matrix::operator=(const Matrix & A)
	{
		for (size_t i = 0; i<rows; i++)
			delete[] vals[i];
		delete[] vals; Matrix();

		rows = A.rows;
		cols = A.cols;
		vals = new double*[rows];
		for (size_t i = 0; i < rows; i++)
			vals[i] = new double[cols];

	}

	Matrix Matrix::operator*(const Matrix & A)
	{
		Matrix B = Matrix(rows, A.cols);
		for (size_t i = 0; i < rows; i++) {
			for (size_t k = 0; k < cols; k++) {
				// double tmp = vals[i][k];
				for (size_t j = 0; j < cols; j++) {
					B.vals[i][j] = vals[i][k] * A.vals[k][j];
				}
			}
		}
		return B;
	}

	Vector Matrix::operator*(const Vector & v)
	{
		Vector x = Vector(rows, 0);
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				x.vals[i] += vals[i][j] * v.vals[j];
			}
		}
		return x;
	}

	Matrix Matrix::operator*(const double & a)
	{
		Matrix B = Matrix(rows, cols);
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				B.vals[i][j] = vals[i][j] * a;
			}
		}
		return B;
	}

	void Matrix::mult(const double & a)
	{
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				vals[i][j] *= a;
			}
		}
	}

	Matrix Matrix::operator+(const Matrix & A)
	{
		Matrix B = Matrix(A.rows, A.cols);
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				B.vals[i][j] = vals[i][j] + A.vals[i][j];
			}
		}
		return B;
	}


	Matrix Matrix::operator-(const Matrix & A)
	{
		Matrix B = Matrix(A.rows, A.cols);
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				B.vals[i][j] = vals[i][j] + A.vals[i][j];
			}
		}
		return B;
	}

	Vector Matrix::operator/(const Vector & b)
	{
		// b.n must equal rows
		double* x = new double[rows];

		double maxA, *ptr;
		size_t imax, i, j, k;

		size_t * p = new size_t[rows];
		double** a = new double*[rows];
		for (i = 0; i < rows; i++) {
			p[i] = i;
			a[i] = new double[cols];
			std::memcpy(a[i], vals[i], cols * 8);
		}

		for (i = 0; i < rows; ++i) {
			maxA = 0.0;
			imax = i;

			for (k = i; k < rows; k++)
				if (abs(a[k][i]) > maxA) {
					maxA = abs(a[k][i]);
					imax = k;
				}

			if (maxA < 1e-64) return false;

			if (imax != i) {
				//pivoting P
				j = p[i];
				p[i] = p[imax];
				p[imax] = j;

				//pivoting rows of A
				ptr = a[i];
				a[i] = a[imax];
				a[imax] = ptr;
			}

			for (j = i + 1; j < rows; j++) {
				a[j][i] /= a[i][i];
				for (k = i + 1; k < rows; k++)
					a[j][k] -= a[j][i] * a[i][k];
			}
		}

		for (i = 0; i < rows; i++) {
			x[i] = b.vals[p[i]];

			for (k = 0; k < i; k++)
				x[i] -= a[i][k] * x[k];
		}

		for (i = rows - 1; i >= 0; i--) {
			for (k = i + 1; k < rows; k++)
				x[i] -= a[i][k] * x[k];

			x[i] /= a[i][i];
		}
		delete[] p;
		for (i = 0; i<rows; i++)
			delete[] a[i];
		delete[] a;

		return Vector(x, rows);
	}

	// Inverse
	Matrix Matrix::operator!()
	{
		Matrix B = Matrix(rows, cols);
		if (rows == 2) {
			double d = 1 / (vals[0][0] * vals[1][1] - vals[1][0] * vals[0][1]);
			B.vals[0][0] = vals[1][1] * d;
			B.vals[0][1] = -vals[1][0] * d;
			B.vals[1][0] = -vals[0][1] * d;
			B.vals[1][1] = vals[0][0] * d;
			return B;
		}
		else {
			Matrix buffer = Matrix(rows - 1, cols - 1);
			double det = 0;
			double sign = 1;
			for (size_t i = 0; i < rows; i++) {
				for (size_t j = 0; j < cols; j++) {
					if (i + j % 2 == 0) {
						B.vals[i][j] = minorDet(i, j, &buffer);
					}
					else {
						B.vals[i][j] = -minorDet(i, j, &buffer);
					}
					if (i == 0) {
						det += B.vals[0][j] * vals[0][j];
						sign *= -1;
					}
				}
			}
			B.mult(1 / det);
			return +B;
		}

	}

	Matrix Matrix::minor(size_t r, size_t c) {
		Matrix B = Matrix(rows - 1, cols - 1);
		size_t ii = 0;
		size_t jj = 0;
		for (size_t i = 0; i < rows; i++) {
			if (i != r) {
				for (size_t j = 0; j < cols; j++) {
					if (j != c) {
						B.vals[ii][jj] = vals[i][j];
						jj++;
					}

				}
				ii++;
			}
		}
		return B;
	}

	void Matrix::minor(size_t r, size_t c, Matrix * B) {
		int ii = 0;
		int jj = 0;
		for (size_t i = 0; i < rows; i++) {
			if (i != r) {
				for (size_t j = 0; j < cols; j++) {
					if (j != c) {
						B->vals[ii][jj] = vals[i][j];
						jj++;
					}

				}
				ii++;
			}
		}
	}

	double Matrix::minorDet(size_t r, size_t c, Matrix * B) {
		size_t ii = 0;
		size_t jj = 0;
		for (size_t i = 0; i < rows; i++) {
			if (i != r) {
				for (size_t j = 0; j < cols; j++) {
					if (j != c) {
						B->vals[ii][jj] = vals[i][j];
						jj++;
					}

				}
				ii++;
			}
		}
		return B->det();
	}

	// transpose
	Matrix Matrix::operator~()
	{
		Matrix B = Matrix(cols, rows);
		for (size_t i = 0; i < rows; i++)
			for (size_t j = 0; j < cols; j++)
				B.vals[j][i] = vals[i][j];
		return B;
	}

	// transpose
	Matrix Matrix::operator+()
	{
		double tmp;
		for (size_t i = 0; i<rows; i++)
		{
			for (size_t j = i + 1; j<cols; j++)
			{
				tmp = vals[i][j];
				vals[i][j] = vals[j][i];
				vals[j][i] = tmp;
			}
		}
		return *this;
	}

	Matrix Matrix::operator*=(const Matrix & A)
	{
		Matrix B = Matrix(rows, cols);
		for (size_t i = 0; i < rows; i++) {
			for (size_t j = 0; j < cols; j++) {
				B.vals[i][j] = A.vals[i][j] * vals[i][j];
			}
		}
		return B;
	}

	Matrix Matrix::mult(const Matrix & A)
	{
		size_t i, j, k, ii, kk;
		size_t ib = 20;
		size_t kb = 20;
		double acc00, acc01, acc10, acc11;
		double** C = new double*[rows];
		for (i = 0; i < rows; i++)
			C[i] = new double[A.cols];
		for (ii = 0; ii < rows; ii += ib)
		{
			for (kk = 0; kk < rows; kk += kb)
			{
				for (j = 0; j < rows; j += 2)
				{
					for (i = ii; i < ii + ib; i += 2)
					{
						if (kk == 0)
							acc00 = acc01 = acc10 = acc11 = 0.0;
						else
						{
							acc00 = C[i + 0][j + 0];
							acc01 = C[i + 0][j + 1];
							acc10 = C[i + 1][j + 0];
							acc11 = C[i + 1][j + 1];
						}
						for (k = kk; k < kk + kb; k++)
						{
							acc00 += vals[k][j + 0] * A.vals[i + 0][k];
							acc01 += vals[k][j + 1] * A.vals[i + 0][k];
							acc10 += vals[k][j + 0] * A.vals[i + 1][k];
							acc11 += vals[k][j + 1] * A.vals[i + 1][k];
						}
						C[i + 0][j + 0] = acc00;
						C[i + 0][j + 1] = acc01;
						C[i + 1][j + 0] = acc10;
						C[i + 1][j + 1] = acc11;
					}
				}
			}
		}
		return Matrix(C, rows, A.cols);
	}

	double Matrix::det()
	{
		// must be square
		if (rows == 2) {
			return vals[1][1] * vals[0][0] - vals[1][0] * vals[0][1];
		}
		else {
			double sum = 0;
			double sign = 1;
			bool * j1 = new bool[cols];
			for (size_t i = 0; i < cols; i++)
				j1[i] = true;

			for (size_t j = 0; j < cols; j++) {
				bool * j2 = new bool[cols];
				std::memcpy(j2, j1, sizeof(j2));
				j2[j] = false;
				sum += sign*vals[0][j] * det(1, j2);
				sign *= -1;
			}
			return sum;
		}
	}

	double Matrix::det(size_t i, bool * j1) {

		if (rows - i > 2) {
			double sum = 0;
			double sign = 1;
			for (size_t j = 0; j < cols; j++) {
				if (j1[j]) {
					bool * j2 = new bool[cols];
					std::memcpy(j2, j1, sizeof(j2));
					j2[j] = false;
					sum += sign*vals[i][j] * det(i++, j2);
					sign *= -1;
				}
			}
			return sum;
		}
		else {
			double a, b;
			size_t j = 0;
			while (!j1[j]) { j++; }
			a = vals[i][j];
			b = vals[i + 1][j];
			while (!j1[++j]) {}
			a = a*vals[i + 1][j] - b*vals[i][j];
			return a;
		}
	}










}