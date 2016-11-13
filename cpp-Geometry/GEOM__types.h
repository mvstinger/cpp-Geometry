/*
 * GEOM__types.h
 *
 *  Created on: Nov 13, 2016
 *      Author: mvstinger
 */

#ifndef GEOM__TYPES_H_
#define GEOM__TYPES_H_



namespace Geometry {



const double GEOM__MAX_ERROR = 1e-6;


const int GEOM__INDETERMINATE =		-100;
const int GEOM__NO_INTERSECTION =	+100;
const int GEOM__INTERSECTION =		+101;



//	Forward declarations
struct Coordinate;
struct CartesianCoordinate;
struct CylindricalCoordinate;
struct SphericalCoordinate;
struct Direction;

class Point;
class Line;
class Segment;
class Vector;



struct Coordinate {};

struct CartesianCoordinate : Coordinate {
	double x;
	double y;
	double z;
	CartesianCoordinate(void);
	CartesianCoordinate(const double, const double, const double);
	CartesianCoordinate(const CartesianCoordinate&);
	CartesianCoordinate(const CylindricalCoordinate&);
	CartesianCoordinate(const SphericalCoordinate&);
};

struct CylindricalCoordinate : Coordinate {
	double rho;
	double phi;
	double z;
	CylindricalCoordinate(void);
	CylindricalCoordinate(const double, const double, const double);
	CylindricalCoordinate(const CartesianCoordinate&);
	CylindricalCoordinate(const CylindricalCoordinate&);
	CylindricalCoordinate(const SphericalCoordinate&);
};

struct SphericalCoordinate : Coordinate {
	double r;
	double theta;
	double phi;
	SphericalCoordinate(void);
	SphericalCoordinate(const double, const double, const double);
	SphericalCoordinate(const CartesianCoordinate&);
	SphericalCoordinate(const CylindricalCoordinate&);
	SphericalCoordinate(const SphericalCoordinate&);
};



struct Direction {
	double theta;
	double phi;

	Direction(void);
	Direction(const double, const double);
	Direction(const Direction&);
};



class Point {
public:
	Point(void);
	Point(const double, const double, const double);
	Point(const CartesianCoordinate&);
	Point(const CylindricalCoordinate&);
	Point(const SphericalCoordinate&);

	double x(void) const;
	double y(void) const;
	double z(void) const;
	double cyl_rho(void) const;
	double cyl_phi(void) const;
	double cyl_z(void) const;
	double sph_r(void) const;
	double sph_theta(void) const;
	double sph_phi(void) const;

	Vector as_Vector(void) const;

private:
	CartesianCoordinate cart_;
	CylindricalCoordinate cyl_;
	SphericalCoordinate sph_;
};



class Line {
public:
	Line(void);
	Line(const Point, const Point);

	Point point_1;
	Point point_2;

	//TODO: define line-line intersection
	int intersection(const Line, Point&) const;
	//TODO: define line-segment intersection
//	int intersect(const Segment);
};


//TODO: de-inherit from Line
class Segment : public Line {
public:
	Segment(void);
	Segment(const Point, const Point);
	Segment(const Segment&);

	double length(void) const;

	Vector as_Vector(void) const;
	explicit operator Vector(void) const;
};


//TODO: de-inherit from Line
class Vector : public Line {
public:
	Vector(void);
	Vector(const Point);
	Vector(const Direction, const double);
	Vector(const Vector&);

	Segment as_Segment(void) const;
	explicit operator Segment(void) const;

	double length(void) const;
	//TODO: Fix norm to be perpendicular
	Vector norm(void) const;

	double x(void) const;
	double y(void) const;
	double z(void) const;

	Vector operator+(const Vector&) const;
	Vector operator-(const Vector&) const;
	double operator*(const Vector&) const;
	Vector operator*(const double) const;
	Vector operator/(const double) const;
	Vector cross(const Vector&) const;
};
inline Vector operator*(const double, const Vector);



const Point ORIGIN = Point(0.0, 0.0, 0.0);



};


#endif /* GEOM__TYPES_H_ */
