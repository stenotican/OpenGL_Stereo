#pragma once

#include "IntPair.h"
#include "cv3.h"

#ifndef _NU_OPTIMIERUNG_H_
  #define _NU_OPTIMIERUNG_H_

// 

//
class nu_Functor1
{
public:


public:
	nu_Functor1(void){}
	virtual ~nu_Functor1(void){}
public:
	virtual double w(double t){t;return 0.0;} 
	virtual double operator()(double t){t;return 0.0;} 
};
//
class nu_Functor2
{
public:
	nu_Functor2(void){}
	virtual ~nu_Functor2(void){}
public:
	virtual double operator()(double t1, double t2){t1;t2;return 0.0;} 
};
//

//
//
class nu_FunctorN
{
public:
	nu_FunctorN(void){}
	virtual ~nu_FunctorN(void){}
public:
	virtual double operator()(CDoubleVector& vec){vec;return 0.0;} 
};
//
class nu_FunctorN_WurmlochTW :public nu_FunctorN
{
public:
	double t1, w1, t2, w2;
	int iNZwischen; //vec.size() = 2*iNZwischen
	virtual double nu_FunctorN_WurmlochTW::operator()(CDoubleVector& vec);

};
//
class nu_FunctorN_Wurmloch :public nu_FunctorN
{
public:
	double t1, w1, t2, w2;
	int N; //vec.size() = N;
	std::vector < cv3> v_cv3;
public:
	double getT(int i);//f(-1) = t1, f(N) = t2; z=t/sqrt(1+t*t); r=sqrt(1+t*t)
	double getW(int i, CDoubleVector& vec);
	void initVec(int iArt, CDoubleVector& vec);

	virtual double operator()(CDoubleVector& vec);

};
//
//
class nu_FunctorN_Wurmloch_Polygon :public nu_FunctorN
{
public:
	double t1, w1, t2, w2;
	int N; //vec.size() = N; = Polygongrad - 1
	int NLines; //Anzahl der Segmente

public:
	double polygon(double t, CDoubleVector& vec, double w1, double w2);
	virtual double operator()(CDoubleVector& vec);

};
//
//

//
class nu_FunctorN_MinKugel :public nu_FunctorN
{
public:
	int m_iOptiZiel;
	int NN; //Anzahl kugeln
	std::vector < cv3> v_cv3;
public:
	nu_FunctorN_MinKugel(){ m_iOptiZiel = 0; }
	void initVec(int iArt, CDoubleVector& vec);
	void reset_vec_cv3(CDoubleVector& vec);
	virtual double operator()(CDoubleVector& vec);
	double mindist();
};
//
////////////////////////////////////////////////////////////////////////
//
//
class nu_Optimierung
{
public:
	int NMAX;
	int nFunktionsAufrufe;

	double phi;
  double resphi;

public:
	nu_Optimierung(void);
	virtual ~nu_Optimierung(void);
public:
	void nrerror(System::String^ s){}
public:
	bool   GoldenerSchnitt_Start(nu_Functor1& f, double& x1, double& x3, int N=10);
	double GoldenerSchnittR(nu_Functor1& f, double& x1, double& x2, double& x3, double eps);
	double GoldenerSchnitt(nu_Functor1& f, double& x1, double& x3, double eps);
//
//
void make_amoeba_StartMatrix(CDoubleMatrix & p, 
														 CDoubleVector &p0, int NDim, double lambda);
//
bool nu_amoeba(int NDim,double ftol,double lambda,
							 CDoubleVector start,
							 nu_FunctorN& f,
							 CDoubleVector& erg);

//
bool nu_amoeba_xmal(int NDim,double ftol,double lambda,
							 CDoubleVector start,
							 nu_FunctorN& f,
							 CDoubleVector& erg);

//
void nu_amoeba(CDoubleMatrix &p, CDoubleVector &y,
							 const double ftol, nu_FunctorN& funk,
	             int &nfunk);

//
double  nu_amotry(CDoubleMatrix &p, 
									CDoubleVector &y, 
									CDoubleVector &psum,
									nu_FunctorN& funk,
	                const int ihi, const double fac);
};
//

//

#endif
