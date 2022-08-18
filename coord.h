//coord.h
#ifndef _COORD_H_INCLUDED_
#define _COORD_H_INCLUDED_
//======================
//forward dependencies

//======================
//include dependencies
#include <iostream>
using namespace std;
//======================

class Coord {

	protected:
		double x;
		double y;
	
	public:
		// Constructors
		Coord(); // default => sets it to (0,0)
		Coord(double x, double y);
		Coord(const Coord& c);
		// Getter Functions
		double get_X() const;
		double get_Y() const;
		// Overloading Operators
		void operator=(const Coord& c);
		Coord operator-(const Coord& c) const;
		Coord operator+(const Coord& c) const;
		void operator+=(const Coord& c);
		void operator-=(const Coord& c);
		Coord operator/(const double d) const;
		Coord operator*(const double d) const;
		Coord distribute(const Coord& c) const;
		bool operator==(const Coord& c); 
		bool operator!=(const Coord&c);
		// Higher level math tasks
		double dot(const Coord& c) const;
		double cross(const Coord& c) const;
		Coord projectOnto(const Coord& c) const;
		Coord perpVector() const;
		double length() const;
		// Display Functions
		friend ostream& operator<<(ostream& os, const Coord& c);


};


#endif
