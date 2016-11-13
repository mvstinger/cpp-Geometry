/*
 * GEOM__types.cc
 *
 *  Created on: Nov 13, 2016
 *      Author: mvstinger
 */

#include <cmath>
#include "../GEOM__types.h"



namespace Geometry {



CartesianCoordinate::CartesianCoordinate(void) :
		x(0.0),
		y(0.0),
		z(0.0) {};

CartesianCoordinate::CartesianCoordinate(const double x, const double y, const double z) :
		x(x),
		y(y),
		z(z) {};

CartesianCoordinate::CartesianCoordinate(const CartesianCoordinate& cart) :
		x(cart.x),
		y(cart.y),
		z(cart.z) {};

CartesianCoordinate::CartesianCoordinate(const CylindricalCoordinate& cyl) :
		x(cyl.rho * cos(cyl.phi)),
		y(cyl.rho * sin(cyl.phi)),
		z(cyl.z) {};

CartesianCoordinate::CartesianCoordinate(const SphericalCoordinate& sph) :
		x(sph.r * sin(sph.theta) * cos(sph.phi)),
		y(sph.r * sin(sph.theta) * sin(sph.phi)),
		z(sph.r * cos(sph.theta)) {};



CylindricalCoordinate::CylindricalCoordinate(void) :
		rho(0.0),
		phi(0.0),
		z(0.0) {};

CylindricalCoordinate::CylindricalCoordinate(const double rho, const double phi, const double z) :
		rho(rho),
		phi(phi),
		z(z) {};

CylindricalCoordinate::CylindricalCoordinate(const CartesianCoordinate& cart) :
		rho(sqrt(cart.x*cart.x + cart.y*cart.y + cart.z*cart.z)),
		phi(atan2(cart.y, cart.x)),
		z(cart.z) {};

CylindricalCoordinate::CylindricalCoordinate(const CylindricalCoordinate& cyl) :
		rho(cyl.rho),
		phi(cyl.phi),
		z(cyl.z) {};

CylindricalCoordinate::CylindricalCoordinate(const SphericalCoordinate& sph) :
		rho(sph.r * sin(sph.theta)),
		phi(sph.phi),
		z(sph.r * cos(sph.theta)) {};



SphericalCoordinate::SphericalCoordinate(void) :
		r(0.0),
		theta(0.0),
		phi(0.0) {};

SphericalCoordinate::SphericalCoordinate(const double r, const double theta, const double phi) :
		r(r),
		theta(theta),
		phi(phi) {};

SphericalCoordinate::SphericalCoordinate(const CartesianCoordinate& cart) :
		r(sqrt(cart.x*cart.x + cart.y*cart.y + cart.z*cart.z)),
		theta(acos(cart.z / sqrt(cart.x*cart.x + cart.y*cart.y + cart.z*cart.z))),
		phi(atan2(cart.y, cart.x)) {};

SphericalCoordinate::SphericalCoordinate(const CylindricalCoordinate& cyl) :
		r(sqrt(cyl.rho*cyl.rho + cyl.z*cyl.z)),
		theta(atan(cyl.rho / cyl.z)),
		phi(cyl.phi) {};

SphericalCoordinate::SphericalCoordinate(const SphericalCoordinate& sph) :
		r(sph.r),
		theta(sph.theta),
		phi(sph.phi) {};



Direction::Direction(void) :
		theta(0.0),
		phi(0.0) {};

Direction::Direction(const double theta, const double phi) :
		theta(theta),
		phi(phi) {};

Direction::Direction(const Direction& dir):
		theta(dir.theta),
		phi(dir.phi) {};



Point::Point(void) :
		cart_(CartesianCoordinate()),
		cyl_(CylindricalCoordinate()),
		sph_(SphericalCoordinate()) {};

Point::Point(const double x, const double y, const double z) :
		cart_(CartesianCoordinate(x, y, z)),
		cyl_(CylindricalCoordinate(CartesianCoordinate(x, y, z))),
		sph_(SphericalCoordinate(CartesianCoordinate(x, y, z))) {};

Point::Point(const CartesianCoordinate& cart) :
		cart_(CartesianCoordinate(cart)),
		cyl_(CylindricalCoordinate(cart)),
		sph_(SphericalCoordinate(cart)) {};

Point::Point(const CylindricalCoordinate& cyl) :
				cart_(CartesianCoordinate(cyl)),
				cyl_(CylindricalCoordinate(cyl)),
				sph_(SphericalCoordinate(cyl)) {};;

Point::Point(const SphericalCoordinate& sph) :
				cart_(CartesianCoordinate(sph)),
				cyl_(CylindricalCoordinate(sph)),
				sph_(SphericalCoordinate(sph)) {};;

double Point::x(void) const { return this->cart_.x; }

double Point::y(void) const { return this->cart_.y; }

double Point::z(void) const { return this->cart_.z; }

double Point::cyl_rho(void) const { return this->cyl_.rho; }

double Point::cyl_phi(void) const { return this->cyl_.phi; }

double Point::cyl_z(void) const { return this->cyl_.z; }

double Point::sph_r(void) const { return this->sph_.r; }

double Point::sph_theta(void) const { return this->sph_.theta; }

double Point::sph_phi(void) const { return this->sph_.phi; }

Vector Point::as_Vector(void) const {
	return Vector( Point(*this) );
}



Line::Line(void) :
		point_1(ORIGIN),
		point_2(ORIGIN) {};

Line::Line(const Point pt1, const Point pt2) :
		point_1(Point(pt1)),
		point_2(Point(pt2)) {};

int Line::intersection(const Line that, Point& intersection) const {
	int intersection_state = GEOM__INDETERMINATE;

	double a1 = this->point_2.x() - this->point_1.x();
	double a2 = that.point_2.x() - that.point_1.x();
	double b1 = this->point_2.y() - this->point_1.y();
	double b2 = that.point_2.y() - that.point_1.y();
	double c1 = this->point_2.z() - this->point_1.z();
	double c2 = that.point_2.z() - that.point_1.z();
	double x1 = this->point_1.x();
	double x2 = that.point_1.x();
	double y1 = this->point_1.y();
	double y2 = that.point_1.y();
	double z1 = this->point_1.z();
	double z2 = that.point_1.z();

	try {
		double t2 = (
				(y1 - y2) / b2 + (b1 / b2) * (x2 - x1) / a1
		) / (
				1 - (b1 / b2) * (a2 / a1) );
		double t1 = (1 / a1) * (x2 - x1) + (a2 / a1) * t2;

		double x_diff = (x1 + a1 * t1) - (x2 + a2 * t2);
		double y_diff = (y1 + b1 * t1) - (y2 + b2 * t2);
		double z_diff = (z1 + c1 * t1) - (z2 + c2 * t2);

		double x_eps = GEOM__MAX_ERROR * abs(x2 - x1);
		double y_eps = GEOM__MAX_ERROR * abs(y2 - y1);
		double z_eps = GEOM__MAX_ERROR * abs(z2 - z1);

		if(x_diff<=x_eps && y_diff<=y_eps && z_diff<=z_eps) {
			intersection_state = GEOM__INTERSECTION;
			intersection = Point(x1 + a1 * t1, y1 + b1 * t1, z1 + c1 * t1);
		} else {
			intersection_state = GEOM__NO_INTERSECTION;
		}
	} catch (...) {
		intersection_state = GEOM__INDETERMINATE;
	}

	return intersection_state;
};


Segment::Segment(void) :
		Line(ORIGIN, ORIGIN) {};

Segment::Segment(const Point pt1, const Point pt2) :
		Line(pt1, pt2) {};

Segment::Segment(const Segment& seg) :
		Line(seg.point_1, seg.point_2) {};

double Segment::length(void) const {
	return sqrt(
			pow(this->point_2.x() - this->point_1.x(), 2) +
			pow(this->point_2.y() - this->point_1.y(), 2) +
			pow(this->point_2.z() - this->point_1.z(), 2) );
};

Vector Segment::as_Vector(void) const {
	return Vector( Point(
			this->point_2.x() - this->point_1.x(),
			this->point_2.y() - this->point_1.y(),
			this->point_2.z() - this->point_1.z() ) );
};

Segment::operator Vector(void) const {
	return this->as_Vector();
};



Vector::Vector(void) :
		Line(ORIGIN, ORIGIN) {};

Vector::Vector(const Point pt2) :
		Line(ORIGIN, ORIGIN) {
	this->point_2 = Point(
				this->point_2.x() - this->point_1.x(),
				this->point_2.y() - this->point_1.y(),
				this->point_2.z() - this->point_1.z() );
};

Vector::Vector(const Direction dir, const double len) :
		Line(ORIGIN, Point(SphericalCoordinate(len, dir.theta, dir.phi))) {};

Vector::Vector(const Vector& vec) :
		Line(vec.point_1, vec.point_2) {};

Segment Vector::as_Segment(void) const {
	return Segment(ORIGIN, Point(this->point_2));
};

Vector::operator Segment() const {
	return this->as_Segment();
};

Vector Vector::norm(void) const {
	return *this / this->length();
};

double Vector::x(void) const {
	return (this->point_2.x() - this->point_1.x());
};

double Vector::y(void) const {
	return (this->point_2.y() - this->point_1.y());
};

double Vector::z(void) const {
	return (this->point_2.z() - this->point_1.z());
};

Vector Vector::operator+(const Vector& that) const {
	return Vector( Point(
			this->x() + that.x(),
			this->y() + that.y(),
			this->z() + that.z() ) );
};

Vector Vector::operator-(const Vector& that) const {
	return Vector( Point(
			this->x() - that.x(),
			this->y() - that.y(),
			this->z() - that.z() ) );
};

double Vector::operator*(const Vector& that) const {
	return (
			this->x() * that.x() +
			this->y() * that.y() +
			this->z() * that.z() );
};

Vector Vector::operator*(const double scalar) const {
	return Vector( Point(
			scalar * this->x(),
			scalar * this->y(),
			scalar * this->z() ) );
};

Vector Vector::operator/(const double scalar) const {
	return Vector( Point(
			this->x() / scalar,
			this->y() / scalar,
			this->z() / scalar ) );
};

Vector Vector::cross(const Vector& that) const {
	return Vector( Point(
			this->y() * that.z() - this->z() * that.y(),
			this->z() * that.x() - this->x() * that.z(),
			this->x() * that.y() - this->y() * that.x() ) );
};

inline Vector operator*(const double scalar, const Vector that) {
	return that * scalar;
};










}



