#include "stdafx.h"
#include <math.h>

#include "nu_optimierung.h"
#include "OpenGL_Base.h"
#include "cv3.h"


using namespace OpenGLForm;
//
double t1, w1, t2, w2;
int iNZwischen;


double nu_FunctorN_WurmlochTW::operator()(CDoubleVector& vec) //dim=2*iZwischen
{
	double ret = 0.0; 

	cv3 v1(t1, w1);
	cv3 v2(t2, w2);

	std::vector<cv3> av;
	av.resize(iNZwischen+2);

	av[0] = v1;
	av[iNZwischen] = v2;
	for (int i = 1; i <= iNZwischen; i++)
	{
		cv3 vi(vec[2 * i - 2], vec[2 * i - 1]);
		av[i] = vi;
	}

	double d1, d2;
	for (int i = 1; i < av.size()-1; i++)
	{
  	d1 = av[i - 1].Dist(av[i]);
		d2 = av[i + 1].Dist(av[i]);
		ret += max(d1, d2);
	}
	return ret;
}

//
//
//
double nu_FunctorN_Wurmloch::getT(int i)//f(-1) = t1, f(N) = t2;
{
	double a = (t2 - t1) / (N + 1.0);
	double b = t1 + a;
	double ret = i*a + b;
	return ret;
}
//
double nu_FunctorN_Wurmloch::getW(int i, CDoubleVector& vec)
{
	if (i < 0)
		return w1;
	if (i >= N)
		return w2;
	return vec[i];
}
//
void nu_FunctorN_Wurmloch::initVec(int iArt, CDoubleVector& vec)
{
	double a = (w2 - w1) / (N + 1.0);
	double b = w1 + a;
	cv3 cv(t1, w1);
	v_cv3.push_back(cv);
	for (int i = 0; i < N; i++)
	{
		vec.setat(i, i*a + b);
		cv = cv3(getT(i), getW(i,vec));
		v_cv3.push_back(cv);
	}
	cv = cv3(t2, w2);
	v_cv3.push_back(cv);
}
//
double nu_FunctorN_Wurmloch::operator()(CDoubleVector& vec)
{
	cv3 cv;
	for (int i = -1; i <= N; i++)
	{
		v_cv3[i + 1].set(getT(i), getW(i, vec));
	}

	double ret = 0.0;
	for (int i = 0; i < v_cv3.size() - 1;i++)
	{
		double d = v_cv3[i].Dist(v_cv3[i + 1]);
		ret += d;
	}
	return ret;
}
//
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//	double t1, w1, t2, w2;
//	int N; //vec.size() = N; = Polygongrad - 1
double nu_FunctorN_Wurmloch_Polygon::polygon(double t, CDoubleVector& vec, double w1, double w2)
{
	double a0 = w1;
	double a1 = w2 - w1 - vec.sum(N);

	int i = N - 1;
	double ret = vec(i);
	while (i > 0)
	{
		i = i - 1;
		ret = t*ret + vec(i);
	}
	ret = t*ret + a1;
	ret = t*ret + a0;
	return ret;
}
//
double nu_FunctorN_Wurmloch_Polygon::operator()(CDoubleVector& vec)
{
	double ret = 0.0;
	v3 v1(t1, w1);
	v3 v2;
	for (int i = 1; i < NLines; i++)
	{
		float t = float(i) / float(NLines);
		double r = polygon(t, vec, t1, t2);
		double w = w1 + t*(w2 - w1);
		v2.set_wurm(r, w);
		ret = ret + v1.Dist(v2);
		v1.set_r(v2);
	}
	return ret;
}
//
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
bool Test_vec_ok(std::vector < cv3>& v)
{
	double dmin = 1.0E20;
	for (int i1 = 0; i1 < v.size(); i1++)
	{
		for (int i2 = i1 + 1; i2 < v.size(); i2++)
		{
			double d = v[i1].Dist(v[i2]);
			if (d < dmin)
				dmin = d;
			if (dmin < 0.25)
				dmin = dmin;
		}
	}
	return (dmin >0.25);
}
//
void nu_FunctorN_MinKugel::reset_vec_cv3(CDoubleVector& vec)
{
	v_cv3.resize(NN);
	v_cv3[0].set_polar(0.0, 0.0);
	v_cv3[1].set_polar(vec[0], 0.0);
	for (int i = 2; i < NN; i++)
	{
		v_cv3[i].set_polar(vec[2 * i - 3], vec[2 * i - 2]);
	}
	Test_vec_ok(v_cv3);
}
//
double nu_FunctorN_MinKugel::mindist()
{
	double dmin = 1.0E20;
	for (int i1 = 0; i1 < NN; i1++)
		for (int i2 = i1 + 1; i2 < NN; i2++)
		{
			double d = v_cv3[i1].Dist(v_cv3[i2]);
			if (d < dmin)
			{
				dmin = d;
			}
		}
	return dmin;
}

//
void nu_FunctorN_MinKugel::initVec(int iArt, CDoubleVector& vec)
{

	double pi = 3.1415926535897932384626433832795;
	double pi2 = 2 * pi;

	switch (iArt)
	{
	case 0:
		{
			vec.setat(0, pi);
			int N1 = (NN - 2) / 2;
			int N2 = NN - 2 - N1;

			double dw = 0.6;
			int N3 = 0;
			for (int i = 0; i < N1; i++)
			{
				vec.setat(2 * i + 1, 0.5*pi - dw); //hw
				vec.setat(2 * i + 2, i*pi2 / N1); //hw
				N3 = 2 * i + 2 + 1;
			}
			for (int i = 0; i < N2; i++)
			{
				vec.setat(2 * i + N3, 0.5*pi + dw); //hw
				vec.setat(2 * i + N3 + 1, pi / N2 + i*pi2 / N2); //hw
			}
			break;
		}
	default:
	{
		int iKreise = iArt / 2;
		int iODD = iArt - 2 * iKreise;
		iKreise++;
		double hw_akt = pi / (iKreise+1);
		int iPktaufKreis = (NN - 2 * iODD) / iKreise;
		int iMittelKreis = (iKreise + 1) / 2;
		int ipktGes = (iODD==1?2:1) + iPktaufKreis*iKreise;
		int ipkt_zusätzlich_auf_Mitte = NN - ipktGes;

		int inext = 0;
		if (iODD == 1)
		{
			vec.setat(0, pi);
			inext = 1;
		}
		else
		{
			vec.setat(0, pi - hw_akt*0.75);
			inext = 1;
		}
		
		double step = 0.0;
		for (int ik = 1; ik <= iKreise; ik++)
		{
			int ik_anz = iPktaufKreis;
			if (ik == iMittelKreis) ik_anz += ipkt_zusätzlich_auf_Mitte;
			step += pi / iPktaufKreis;
			double delta = pi2 / ik_anz;
			for (int i = 0; i < ik_anz; i++)
			{
				vec.setat(inext, hw_akt*ik); inext++;
				vec.setat(inext, i*delta + step); inext++;
			}
		}


	}
	break;
	}

	reset_vec_cv3(vec);
}
//
double nu_FunctorN_MinKugel::operator()(CDoubleVector& vec)
{
	reset_vec_cv3(vec);
	double sum = 0.0;

	switch (m_iOptiZiel)
	{
	default:
		for (int i = 0; i < NN; i++)
			for (int j = i + 1; j < NN; j++)
			{
				double d = v_cv3[i].Dist(v_cv3[j]);
				sum += 1.0 / d;
			}
		break;
	case 1:
		sum = 123456.0;
		for (int i = 0; i < NN; i++)
			for (int j = i + 1; j < NN; j++)
			{
				double d = v_cv3[i].Dist(v_cv3[j]);
				if (d < sum)
					sum = d;
			}
		break;
	}
	return sum;
}
//
//
//
nu_Optimierung::nu_Optimierung()
{
	nFunktionsAufrufe = 0;
	NMAX = 50000;
	phi = (1.0 + sqrt(5.0)) / 2.0;
  resphi = 2 - phi;
}

nu_Optimierung::~nu_Optimierung(void)
{}

//
bool nu_Optimierung::GoldenerSchnitt_Start(nu_Functor1& f, double& x1, double& x3, int N)
{
	double t,w, tmin=0.0;
	double wmin = f(x1)+1.0;
	double ds = (x3-x1)/N;
  for(int i=0;i<=N;i++)
	{
		t = x1 + i*ds;
		w = f(t);
		if(w < wmin)
		{
			tmin = t;
			wmin = w;
		}
	}

	double dds = 0.1*ds;
	while(dds < 2.0*ds)
	{
  	x1 = tmin-dds;
	  x3 = tmin+dds;
	  if(f(x1) > wmin && f(x3) > wmin)
		{
		  return true;
		}
		dds += 0.1*ds;
	}

	return false;
}

//

double nu_Optimierung::GoldenerSchnitt(nu_Functor1& f, double& x1, double& x3, double eps)
{
	double w2, w4, x5,w5;
	    if (fabs(x1 - x3) < eps)
			{
        return (x1 + x3) / 2.0;
			}
 
    double x2 = x1 + resphi * (x3 - x1);
		w2 = f(x2);

    // Create a new possible center in the area between x2 and x3, closer to x2
    double x4 = x2 + resphi * (x3 - x2);
		w4 = f(x4);

 while(fabs(x1 - x3) > eps)
 {

	 if (w4 < w2)
		{
				x1 = x2;
				x2 = x4;
				w2 = w4;
				if((x2-x1) < (x3-x2))
				{
				  x4 = x2 + resphi * (x3 - x2);
				  w4 = f(x4);
				}
				else
				{
					x5 = x2 - resphi *(x2-x1);
					w5 = f(x5);
					x4 = x2;
					w4 = w2;
					x2 = x5;
					w2 = w5;
				}
		}
    else
		{
        x3 = x4; 
				if((x2-x1) < (x3-x2))
				{
				  x4 = x2 + resphi * (x3 - x2);
				  w4 = f(x4);
				}
				else
				{
					x5 = x2 - resphi *(x2-x1);
					w5 = f(x5);
					x4 = x2;
					w4 = w2;
					x2 = x5;
					w2 = w5;
				}
		}
 }
 return (x1 + x3) / 2.0;
}
//
//
double nu_Optimierung::GoldenerSchnittR(nu_Functor1& f, double& x1, double& x2, double& x3, double eps)
{
	  if (fabs(x1 - x3) < eps)
		{
      return (x1 + x3) / 2.0;
		}
 
    // Create a new possible center in the area between x2 and x3, closer to x2
    double x4 = x2 + resphi * (x3 - x2);
 
    if (f.w(x4) < f.w(x2))
		{
        return GoldenerSchnittR(f, x2, x4, x3, eps);
		}
    else
		{
        return GoldenerSchnittR(f, x4, x2, x1, eps);
		}
		return (x1 + x3) / 2.0;
}
//
bool nu_Optimierung::nu_amoeba(int NDim,double ftol,double lambda,
							 CDoubleVector start,
							 nu_FunctorN& f,
							 CDoubleVector& erg)
{
	if(nu_amoeba_xmal(NDim, ftol, lambda, start,f, erg))
	{
		start = erg;
		if(nu_amoeba_xmal(NDim, ftol, lambda, start,f, erg))
		{
			lambda = 0.1*lambda;
			start = erg;
			if(nu_amoeba_xmal(NDim, ftol, lambda, start,f, erg))
			{
				return true;
			}
		}
	}
	return false;
}

bool nu_Optimierung::nu_amoeba_xmal(int NDim,double ftol,double lambda,
							 CDoubleVector start,
							 nu_FunctorN& f,
							 CDoubleVector& erg)
{
  CDoubleMatrix m;
	make_amoeba_StartMatrix(m, start, NDim, lambda);
	CDoubleVector y(NDim+1);
	CDoubleVector pt(NDim);
	for(int i=0;i<=NDim; i++)
	{
		m.get_row(i,pt, NDim);
		double dff = f(pt);
		y.setat(i, dff);
	}

	
	nu_amoeba(m, y, ftol, f, nFunktionsAufrufe);
	double test;
	if (nFunktionsAufrufe < NMAX)
	{
		erg = CDoubleVector(NDim);
		for(int i=0;i<=NDim; i++)
		{
			m.get_row(i,pt, NDim);
			test = f(pt);
			erg.add(NDim, pt);
		}
		for(int i=0;i<NDim;i++)
		{
			double d = erg.getat(i);
			d = d /(NDim+1);
			erg.setat(i,d);
		}
		test = f(erg);
    return true; 
	}
  return false;
}

//
void nu_Optimierung::make_amoeba_StartMatrix(CDoubleMatrix & m, 
														 CDoubleVector &p0, 
														 int NDim, double lambda)
{
	m.set_rc(NDim+1, NDim);
	m.set_row(0, p0, NDim);
	CDoubleVector p;
	for(int i=1;i<=NDim; i++)
	{
		p = p0;
	  p.add(i-1, lambda);
		m.set_row(i, p, NDim);
	}
}
//
	 void get_psum(CDoubleMatrix &p, CDoubleVector &psum)
	{
		int i,j;
		double sum;

		int mpts=p.nrows();
		int ndim=p.ncols();
		for (j=0;j<ndim;j++) 
		{
			for (sum=0.0,i=0;i<mpts;i++)
			{
				sum += p(i,j);
			}
			psum.setat(j,sum);
		}
	}


//
void nu_Optimierung::nu_amoeba(CDoubleMatrix &p, CDoubleVector &y,
															 const double ftol, nu_FunctorN& funk,
	                             int &nfunk)
{
	const double TINY=1.0e-12;
	int i,ihi,ilo,inhi,j;
	double rtol,ysave,ytry;

	int mpts=p.nrows();
	int ndim=p.ncols();
	CDoubleVector psum(ndim);
	nfunk=0;
	get_psum(p,psum);
	for (;;) {
		ilo=0;
		//ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);


		if(y(0)>y(1))
		{
			inhi = 1;
			ihi  = 0;
		}
		else
		{
			inhi = 0;
			ihi  = 1;
		}

		for (i=0;i<mpts;i++) {
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi) inhi=i;
		}
		rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
		if (rtol < ftol) 
		{
			//SWAP(y[0],y[ilo]);
			y.swap(0,ilo);
			for (i=0;i<ndim;i++) p.swap(0,i,ilo,i);//SWAP(p[0][i],p[ilo][i]);
			break;
		}
		if (nfunk >= NMAX)
    {
      nrerror("NMAX exceeded");
      return;
    }
		nfunk += 2;
		ytry=nu_amotry(p,y,psum,funk,ihi,-1.0);
		if (ytry <= y[ilo])
			ytry=nu_amotry(p,y,psum,funk,ihi,2.0);
		else if (ytry >= y[inhi]) {
			ysave=y[ihi];
			ytry=nu_amotry(p,y,psum,funk,ihi,0.5);
			if (ytry >= ysave) {
				for (i=0;i<mpts;i++) {
					if (i != ilo) {
						for (j=0;j<ndim;j++)
						{
							double w=0.5*(p(i,j)+p(ilo,j));
							p.setat(i,j,w);
							psum.setat(j,w);
						}
						y[i]=funk(psum);
					}
				}
				nfunk += ndim;
				get_psum(p,psum);
			}
		} else --nfunk;
	}
}


double  nu_Optimierung::nu_amotry(CDoubleMatrix &p, 
																	CDoubleVector &y, 
																	CDoubleVector &psum, nu_FunctorN& funk,
	const int ihi, const double fac)
{
	int j;
	double fac1,fac2,ytry;

	int ndim=p.ncols();
	CDoubleVector ptry(ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=0;j<ndim;j++)
		ptry.setat(j,psum(j)*fac1-p(ihi,j)*fac2);
	ytry=funk(ptry);
	if (ytry < y(ihi))
	{
		y.setat(ihi,ytry);
		for (j=0;j<ndim;j++) {
			psum.add(j,ptry(j)-p(ihi,j));
			p.setat(ihi,j,ptry(j));
		}
	}
	return ytry;
}
