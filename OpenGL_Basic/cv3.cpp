#include "stdafx.h"
#include "cv3.h"

#include "math.h"

double cv3::eps = 1.0E-6;

cv3::cv3(double r, double w)
{
	double t2 = 1.0 + r*r;
	double t1 = sqrt(t2);
	x = t1*cos(w);
	y = t1*sin(w);
	z = r / t1;
}


cv3::cv3(int iArt, double r, double w) //xyz-Normale
{
	switch (iArt)
	{
	case 0:
	default:
	{
		double t2 = 1 + r*r;
		x = -cos(w)*(t2);
		y = -sin(w)*(t2);
		z = r;
		double n = Norm();
		x = x / n;
		y = y / n;
		z = z / n;
	}
  break;
	case 1://df/dr;
	{
		double t1 = sqrt(1.0 + r*r);
		double t2 = r / t1;
		x = t2*cos(w);
	  y = t2*sin(w);
		z = 1.0 / (t1*t1*t1);
	}
	break;
	case 2://df/dw;
	{
		double t1 = sqrt(1.0 + r*r);
		x = -t1*sin(w);
		y = t1*cos(w);
		z = 0.0;
	}
	break;
	}
}
//
void cv3::set_null()
{
	x = y = z = r = 0.0;
}
//
void cv3::set(double r, double w) //xyz-Wurm
{
	double t2 = 1 + r*r;
	double t1 = sqrt(t2);
	x = t1*cos(w);
	y = t1*sin(w);
	z = r / t1;
}
//
void cv3::set_polar(double hw, double gw) // hw=0 ->z=1;
{
	x = sin(hw)*cos(gw);
	y = sin(hw)*sin(gw);
	z = cos(hw);
}
//
void cv3::get_polar(double& hw, double& gw)
{
	cv3 c = *this;
	c.Normieren(1.0);
	hw = acos(z);
	double sin_hw = sin(hw);
	gw = atan2(y, x);

	cv3 test; test.set_polar(hw, gw);
	double d = test.Dist(c);
	if (d > 0.0000000001)
	{
		d = 0.0;
	}
}
//
//
double cv3::Dist(cv3& v)
{
	return sqrt((v.x - x)*(v.x - x) + (v.y - y)*(v.y - y) + (v.z - z)*(v.z - z));
}
//
double cv3::Dist2(cv3& v)
{
	return ((v.x - x)*(v.x - x) + (v.y - y)*(v.y - y) + (v.z - z)*(v.z - z));
}

//
double cv3::Norm() const 
{
	return sqrt(x*x + y*y + z*z);
}
void cv3::Normieren(double laenge)
{
	double nn = Norm();
	if (nn > 0.0)
	{
		x = x*laenge / nn;
		y = y*laenge / nn;
		z = z*laenge / nn;
	}
}
//
double cv3::Max_xyz() const
{
	double dmax = fabs(x);
	if (fabs(y) > dmax) dmax = fabs(y);
	if (fabs(z) > dmax) dmax = fabs(z);
	return dmax;
}



cv3 operator+(const cv3& a, const cv3& b)
{
	cv3 v(a.x + b.x, a.y + b.y, a.z + b.z);
	return v;
}

cv3 operator-(const cv3& a, const cv3& b)
{
	cv3 v(a.x - b.x, a.y - b.y, a.z - b.z);
	return v;
}
//----------------------------------------------------------------



/*
v3^ operator+(const v3% a, const v3% b)
{
v3^v = gcnew v3();
v->x = a.x + b.x;
v->y = a.y + b.y;
v->z = a.z + b.z;
return v;
}
*/

cv3 operator*(const double a, const cv3& b)
{
	cv3 v;
	v.x = a * b.x;
	v.y = a * b.y;
	v.z = a * b.z;
	return v;
}


double operator*(const cv3& a, const cv3& b)
{
	double f = a.x*b.x + a.y*b.y + a.z*b.z;
	return f;
}

cv3 operator%(const cv3& a, const cv3& b)
{
	cv3 v;
	v.x = a.y*b.z - a.z*b.y;
	v.y = a.z*b.x - a.x*b.z;
	v.z = a.x*b.y - a.y*b.x;
	return v;
}
//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
bool solve2(double a11, double a12, double a21, double a22, double b1, double b2, double& x1, double& x2)
{
	double det = a11*a22 - a12*a21;
	if (fabs(det) > 1.0E-8)
	{
		double det2 = a11*b2 - a21*b1;
		double det1 = a22*b1 - a12*b2;
		x1 = det1 / det;
		x2 = det2 / det;
		return true;
	}
	return false;
}
//----------------------------------------------------------------
double verbindungskante(cv3 p1, cv3 p2, cv3 q1, cv3 q2, double& tp, double& ts, cv3& px, cv3& qx)
{
	cv3 p3 = p2 - p1;
	cv3 q3 = q2 - q1;

	double a11 = p3*p3;
	double a21 = q3*p3;
	double a12 = -a21;
	double a22 = q3*q3; a22 = -a22;
	double b1 = q1*p3 - p1*p3;
	double b2 = q1*q3 - p1*q3;
	double x1, x2;
	if (solve2(a11, a12, a21, a22, b1, b2, x1, x2))
	{
		tp = x1;
		ts = x2;
		px = p1 + tp*p3;
		qx = q1 + ts*q3;
		return px.Dist(qx);
	}
	return _DMAX;
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

/*
vec_cv3 m_pkt;
	vec_cv3 m_ri;
	std::vector<double> m_d01;
	*/

	void bezier::Berechnung()
	{
		m_d01.push_back(0.0);
		for (int i = 1; i < m_pkt.size(); i++)
		{
			double d = m_pkt[i - 1].Dist(m_pkt[i]);
			m_d01.push_back(d);
		}
	}
	cv3 bezier::wert(double d01)
	{
		return cv3();
	}
	////////////////////////////////////////////////////////////////////////
	double winkel(const cv3& a, const cv3& b)
	{
		double na = a.Norm();
		double nb = b.Norm();
		if ((na>0.0) && (nb > 0.0))
		{
			double cosw = (a*b) / (na*nb);
			double ret = acos(cosw);
			return ret;
		}
		return -1.0;
	}
	//
	cv3 cv3::drehung_vektor(cv3 v0, cv3 achse, double rad_winkel)
	{
		double vl = v0.Norm();
		v0.Normieren(1.0);
		cv3 vq = achse%v0;
		vq.Normieren(1.0);
		cv3 ret = cos(rad_winkel)*v0 + sin(rad_winkel)*vq;
		ret = vl*ret;
		return ret;
	}
	//
	cv3 cv3::drehung_punkt(cv3 p0, cv3 achse, double rad_winkel, cv3 mittelpunkt)
	{
		cv3 v0 = p0 - mittelpunkt;
		cv3 ret = drehung_vektor(v0, achse, rad_winkel);
		ret = ret + mittelpunkt;
		return ret;
	}
	//
	cv3 operator- (const cv3& a)
	{
		cv3 ret(a);
		ret.x = -ret.x;
		ret.y = -ret.y;
		ret.z = -ret.z;
		return ret;
	}
	////////////////////////////////////////////////////////////////////////
	cgerade::cgerade(cv3 _p0, cv3 _ri)
	{
		p0 = _p0;
		richtung = _ri;
	}
	cv3 cgerade::punkt01(double t) const
	{
		cv3 p = p0 + t*richtung;
		return p;
	}
	bool cgerade::normalform()
	{
		double rr = richtung*richtung;
		if (rr > 0.00000001)
		{
			double t = -(p0*richtung) / rr;
			p0 = p0 + t*richtung;
			richtung.Normieren(1.0);

			return true;
		}
		return false;
	}
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
	cebene::cebene(cv3 p1, cv3 p2, cv3 p3, bool check_ab_positiv)
	{
		cv3 q1 = p2 - p1;
		cv3 q2 = p3 - p1;
		lot = q1%q2;
		lot.Normieren(1.0);
		ab = lot*p1;
		if (fabs(ab) < 0.05)
		{
			ab = lot*p1;
		}
		if (check_ab_positiv && (ab < 0.0))
		{
			lot = -lot;
			ab = -ab;
		}
	}

	//
	cebene::cebene(cv3 _lot, cv3 p0)
	{
		lot = _lot;
		lot.Normieren(1.0);
		ab = lot*p0;
	}
	//
	double cebene::abstand(cv3 p) const 
	{
		return p*lot - ab;
	}
	////////////////////////////////////////////////////////////////////////
	cv3 schnittpunkt(const cebene& e, const cgerade& g, bool& ok)
	{
		cv3 schnitt;
		double d0 = e.abstand(g.punkt01(0.0));
		double d1 = e.abstand(g.punkt01(1.0));
		double dd = d0 - d1;
		if (fabs(dd) > cv3::eps)
		{
			double t = d0 / dd;
			schnitt = g.punkt01(t);
			ok = true;
		}
		ok = false;
		return schnitt;
	}
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
