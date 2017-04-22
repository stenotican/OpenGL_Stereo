#pragma once

#include <vector>


class cv3
{
public:
	static double  eps;
public:
	double  x, y, z, r;

	public:
		cv3(){ x = 0.0; y = 0.0; z = 0.0; r = 0.0; }
		cv3(const cv3& v){ x = v.x; y = v.y; z = v.z; r = v.r; }
		cv3(double a, double b, double c){ x = a; y = b; z = c; r = 0.0; }
		cv3(double a, double b, double c, double rr){ x = a; y = b; z = c; r = rr; }
		cv3(double r, double w); //xyz-Wurm
		cv3(int iArt, double r, double w); //xyz-Normale
		virtual ~cv3(){}

		void set_null();
		void set(double r, double w); //xyz-Wurm
		void set_polar(double hw, double gw); // hw=0 ->z=1;
		void get_polar(double& hw, double& gw); // umwandlung x,y,z ingw,hw;

	public:
		double Dist(cv3& v);
		double Dist2(cv3& v);
		double Norm()const;// { return sqrt(x*x + y*y + z*z); }
		void Normieren(double laenge);
		double Max_xyz() const;

		cv3& operator=(const cv3& v){ x = v.x; y = v.y; z = v.z; return *this; }

		static cv3 drehung_vektor(cv3 v0, cv3 achse, double rad_winkel);
		static cv3 drehung_punkt(cv3 v0, cv3 achse, double rad_winkel, cv3 mittelpunkt);
};

double winkel(const cv3& a, const cv3& b);

cv3 operator- (const cv3& a);

cv3 operator+(const cv3& a, const cv3& b);
cv3 operator-(const cv3& a, const cv3& b);
cv3 operator*(const double a, const cv3& b);

double operator*(const cv3& a, const cv3& b);
cv3 operator%(const cv3& a, const cv3& b);

//geometrie
bool solve2(double a11, double a12, double a21, double a22, double b1, double b2, double& x1, double& x2);
double verbindungskante(cv3 p1, cv3 p2, cv3 q1, cv3 q2, double& tp, double& ts, cv3& p3, cv3&q3);


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

class vec_cv3 :public std::vector < cv3 >
{
public:
	void add(double x, double y, double z){ cv3 v(x, y, z);  push_back(v); }
	void add(double x, double y, double z, double r){ cv3 v(x, y, z,r);  push_back(v); }
};
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
class cgerade
{
public:
	cv3 p0, richtung;
public:
	cgerade(cv3 _p0, cv3 _ri);

public:
	cv3 punkt01(double t) const;
	bool normalform();

};

///////////////////////////////////////////////////////////////////////////////////////
class cebene
{
public:
	cv3 lot;
	double ab;
public:
	cebene(cv3 p1, cv3 p2, cv3 p3, bool check_ab_positiv=false);
	cebene(cv3 _lot, cv3 p0);
	double abstand(cv3 p) const;
};

cv3 schnittpunkt(const cebene& e, const cgerade& g, bool& ok);
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

class bezier
{
public:
	vec_cv3 m_pkt;
	vec_cv3 m_ri;
	std::vector<double> m_d01;
public:
	void Berechnung();
	cv3 wert(double d01);

};