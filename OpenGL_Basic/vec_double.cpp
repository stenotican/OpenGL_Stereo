#include "stdafx.h"
#include "vec_double.h"
#include <vector>
#include <string>
#include <algorithm>

#define Try try{
#define Catch }catch(...){}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void nrerror(System::String^ msg)
{
 
}

vec_int::vec_int(){};

vec_int::vec_int(int N)
{
  for(int i = 0; i<N; i++)
  {
    push_back(0);
  }
}

vec_int::vec_int(const vec_int& objectSrc)
{
  this->operator =(objectSrc);
}

vec_int& vec_int::operator =(const vec_int& objectSrc)
{
  this->std::vector<int>::operator =(objectSrc);
	//(*this).operator =(objectSrc)
  return *this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
vec_double::vec_double(){};

vec_double::vec_double(const double&x)
{
  this->push_back(x);
}

vec_double::vec_double(const double&x, const double&y)
{
  this->push_back(x);
  this->push_back(y);
}

vec_double::vec_double(const double&x, const double&y, const double&z)
{
  this->push_back(x);
  this->push_back(y);
  this->push_back(z);
}

vec_double::vec_double(const double&x, const double&y, const double&z, const double&w)
{
  this->push_back(x);
  this->push_back(y);
  this->push_back(z);
  this->push_back(w);
}


vec_double::vec_double(int N)
{
  for(int i = 0; i<N; i++)
  {
    push_back(0.0);
  }
 
}

vec_double::vec_double(const vec_double& objectSrc)
{
  this->operator =(objectSrc);
}

vec_double& vec_double::operator =(const vec_double& objectSrc)
{
  this->std::vector<double>::operator =(objectSrc);
  //(*this).operator =(objectSrc);
  return *this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double fmin(double a, double b)
{
  if (a < b)
    return a;
  return b;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
vec_double operator+(const vec_double& a, const vec_double& b)                                  // Addition        // freie Funktion
{
	Try
    vec_double c;
    int min_array = fmin(a.size(), b.size());
    for(int i = 0; i < min_array; i++)
    {
      c.push_back(a[i] + b[i]);
    }
    return c;
  Catch;
  return vec_double();
  nrerror("Fehler:vec_double::operator+");
}

vec_double operator-(const vec_double& a, const vec_double& b)                                  // Subtraktion      // freie Funktion
{
	Try
    vec_double c;
    int min_array = fmin(a.size(), b.size());
    for(int i = 0; i < min_array; i++)
    {
      c.push_back(a[i] - b[i]);
    }
    return c;
  Catch;
  return vec_double();
  nrerror("Fehler:vec_double::operator-");
}

double operator*(const vec_double& a, const vec_double& b)                                      //  Skalarprodukt(x*x)+(y*y)+(z*z)...     // freie Funktion
{
	Try
    double skalarprodukt = 0.0;
    int min_array = fmin(a.size(), b.size());
    for(int i = 0; i < min_array; i++)
    {
      skalarprodukt = skalarprodukt +(a[i] * b[i]);
    }
    return skalarprodukt;
  Catch
  nrerror("Fehler:operator*");
  return 0.0;
}

vec_double operator-(const vec_double& a)                                                       // Vorzeichenwechsel    // freie Funktion
{
	Try
    vec_double c;
    for(int i = 0; i < a.isize(); i++)
    {
      c.push_back(-a[i]);
    }
    return c;
  Catch;
  nrerror("Fehler:vec_double::operator-");
  return vec_double();
}

vec_double operator*(const double a, const vec_double& b)                                       // ein double *vektor z.B 3.0*(1/2/3)   // freie Funktion
{
	Try
    vec_double c;
    for(int i = 0; i < b.isize(); i++)
    {
      c.push_back(a*(b[i]));
    }
    return c;
  Catch
  nrerror("Fehler:vec_double::operator*");
  return vec_double();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double vec_double::maximum()const
{
	Try
    double max = 0;
    for(int i = 0; i < isize(); i++)
    {
      if((*this)[i] > max)
      {
        max =(*this)[i];
      }
    }
    return max;
  Catch;
  nrerror("Fehler:vec_double::maximum()");
  return 0.0;
}

double vec_double::minimum()const
{
	Try
    double min = 100000000000;
    for(int i = 0; i < isize(); i++)
    {
      if((*this)[i] < min)
      {
        min =(*this)[i];
      }
    }
    return min;
  Catch
  nrerror("Fehler:vec_double::minimum()");
  return 0.0;
}

double vec_double::skalarprodukt(const vec_double& a, const vec_double& b)                      // zwei Vektoren multiplizieren::(x*x)+(y*y)+(z*z)...
{
	Try
    double skalarprodukt = 0.0;
    int min_array = fmin(a.size(), b.size());
    for(int i = 0; i < min_array; i++)
    {
      skalarprodukt = skalarprodukt +(a[i] * b[i]);
    }
    return skalarprodukt;
  Catch
  nrerror("Fehler:vec_double::skalarprodukt");
  return 0.0;
}

double vec_double::norm()const                                                                   // Wurzel aus einen Vektor(x²+y²+z²...)
{
	Try
    double norm = 0.0;
    for(int i = 0; i < isize(); i++)
    {
      norm = norm +((*this)[i] *(*this)[i]);
    }
    norm = sqrt(norm);
    return norm;
  Catch
  nrerror("Fehler:vec_double::norm()");
  return 0.0;
}

void vec_double::normieren(double länge)                                                        // Vektor auf eine bestimmte Länge setzen
{
	Try
    double länge_alt = norm();
    if(länge > 0)
    {
      if(länge_alt != länge)
      {
        double a = länge / länge_alt;
        (*this)= a*(*this);
      }
      if(länge != norm())
      {
        nrerror("wer");
      }
    }
    else
    {
      nrerror("wer2");

    }
  Catch
  nrerror("Fehler:vec_double::normieren");
}

double vec_double::summe()const
{
	Try
    double summe = 0.0;
    for(int i = 0; i < isize(); i++)
    {
      summe = summe +(*this)[i];
    }
    return summe;
  Catch
  nrerror("Fehler:vec_double::summe()");
  return 0.0;
}

double vec_double::abstand(const vec_double&a)const                                              // Abstand zw Vektoren bestimmen
{
	Try
    double abstand = 0.0;
    vec_double differenz =(*this)- a;
    abstand = differenz.norm();
    return abstand;
  Catch
  nrerror("Fehler:vec_double::abstand");
  return 0.0;
}

double vec_double::winkel(const vec_double& a, const vec_double& b)const                         // Winkel zw Vektoren bestimmen
{
	Try
    double winkel1 = 0.0;
    if(a.norm()*b.norm()!= 0)
    {
      double temp = a*b /(a.norm()*b.norm());
      winkel1 = acos(temp);
    }
    else
    {
      winkel1 = 90.0;
    }
    return winkel1;
  Catch
  nrerror("Fehler:vec_double::winkel");
  return 0.0;
}


double vec_double::getat(int i) const
{
  return (*this)[i];
}


void vec_double::setat(int index, double iWert)
{
  if (index >= 0 && index < this->isize())
  {
    (*this)[index] = iWert;
  }  
}


void vec_double::add(int i, double d_add)                         // addieren
{
  double d = getat(i);
  d += d_add;
  setat(i, d);
}

void vec_double::add_vec(int NDim, const vec_double& v)
{
  for (int i = 0; i<NDim; i++)
  {
    add(i, v.getat(i));
  }
}


void vec_double::swap(int i, int j)
{
  double di = getat(i);
  double dj = getat(j);
  setat(i, dj);
  setat(j, di);
}


mat_double::mat_double()
{
  m_nrows = 0;
  m_ncols = 0;
}

mat_double::mat_double(int Zeilen, int Spalten)                                                 // Zeilen und Spalten der Matrix auf 0.0 setzen/auffüllen
{
  m_nrows = Zeilen;
  m_ncols = Spalten;
  m_Data.resize(Zeilen*Spalten, 0.0);
}

mat_double::mat_double(vec_double v, Zeile_Spalte eZS)
{
  switch(eZS)
  {
    case eZeile:
      m_nrows = 1;
      m_ncols = v.isize();
      m_Data.resize(m_nrows*m_ncols, 0.0);
      set_row(0, v, 0);
      break;
    case eSpalte:
      m_ncols = 1;
      m_nrows = v.isize();
      m_Data.resize(m_nrows*m_ncols, 0.0);
      set_col(0, v, 0);
      break;
  }
  if(m_nrows < 1 ||m_ncols< 1)
  {
    nrerror("Fehler:mat_double::mat_double");
  }
}

void mat_double::init()                                                                         // die Matrix mit Werten füllen
{
	Try
    for(int i = 0; i < m_Data.isize(); i++)
    {
      m_Data[i] = 1/double(i+1);
    }
    return;
  Catch
  nrerror("Fehler:mat_double::init()");
}

double mat_double::dummy = 0.0;


double& mat_double::operator()(int iZ, int iS)                                                  //()Operator
{
	Try
    if(iZ >= 0 && iS >= 0 &&iZ<m_nrows &&iS<m_ncols)                                            // Zugrif auf einen Wert in der Matrix(Zeile, Spalte)
    {
      int index = iZ*m_ncols+ iS;
      return m_Data[index];
    }
    else
    {
      nrerror("Fehler:mat_double::operator()");
    }
  Catch
nrerror("Fehler Catch:mat_double::operator()");
return dummy;
}

double* mat_double::operator[](int iZ)                                                          // [] Operator
{
	Try
    if(iZ >= 0 &&  iZ<m_nrows)                                                                  // Zugrif auf eine [Zeile] oder [Spalte]
    {
      int index = iZ*m_ncols ;
      return &m_Data[index];
    }
    else
    {
      nrerror("Fehler:mat_double::operator[]");
    }
  Catch
nrerror("Fehler Catch:mat_double::operator[]");
return &dummy;
}

void mat_double::swap(int z1, int sp1, int z2, int sp2)
{
  double di = get(z1, sp1);
  double dj = get(z2, sp2);
  set(z1, sp1, dj);
  set(z2, sp2, di);
}





void mat_double::set(int ir, int ic, double wert)                                               // bestimmte Stelle in der Matrix mit einem Wert füllen
{
	Try
    if(ir >= 0 && ic >= 0 && ir < m_nrows && ic < m_ncols)
    {
      int index = ir*m_ncols + ic;
      m_Data[index] = wert;
    }
    else
    {
      nrerror("Fehler:mat_double::set");
    }
    return;
  Catch
nrerror("Fehler Catch:mat_double::set");
}
void mat_double::set_row(int ir, vec_double vec, int /*Ndim*/)                                      // einer Zeile den Vector zuweisen
{
	Try
    if(ir >= 0 && ir < m_nrows && vec.isize()<= m_ncols)
    {
      for(int ic = 0; ic < vec.isize(); ic++)
      {
        set(ir, ic, vec[ic]);
      }
    }
    else
    {
      nrerror("Fehler:mat_double::set_row");
    }
    return;
  Catch
nrerror("Fehler Catch:mat_double::set_row");
}
void mat_double::set_col(int ic, vec_double vec, int /*Ndim*/)                                      // einer Spalte den Vector zuweisen
{
	Try
    if(ic >= 0 && ic < m_ncols && vec.isize()<= m_nrows)
    {
      for(int ir = 0; ir < vec.isize(); ir++)
      {
        set(ir, ic, vec[ir]);
      }
    }
    else
    {
      nrerror("Fehler:mat_double::set_col");
    }
    return;
  Catch
nrerror("Fehler Catch:mat_double::set_col");
}
double& mat_double::get(int ir, int ic)                                                         // einen Wert aus der Matrix auslesen
{
	Try
    if(ir >= 0 && ic >= 0 && ir<m_nrows &&ic<m_ncols)
    {
      int index = ir*m_ncols + ic;
      return m_Data[index];
    }
    else
    {
      nrerror("Fehler:mat_double::get");
    }
  Catch
nrerror("Fehler Catch:mat_double::get");
return dummy;
}
void mat_double::get_row(int ir, vec_double& vec, int /*Ndim*/)                                     // Werte aus einer Zeile der Matrix bekommen
{
	Try
    if(ir >= 0 && ir<m_nrows)
    {
      vec.resize(0);
      for(int ic = 0; ic < m_ncols; ic++)
      {
        vec.push_back(get(ir, ic));
      }
    }
    else
    {
      nrerror("Fehler:mat_double::get_row");
    }
    return;
  Catch
nrerror("Fehler Catch:mat_double::get_row");
}
void mat_double::get_col(int ic, vec_double& vec, int /*Ndim*/)                                     // Werte aus einer Spalte der Matrix bekommen
{
	Try
    if(ic >= 0 && ic<m_nrows)
    {
      vec.resize(0);
      for(int ir = 0; ir < m_nrows; ir++)
      {
        vec.push_back(get(ir, ic));
      }
    }
    else
    {
      nrerror("Fehler:mat_double::get_col");
    }
    return;
    Catch
}

void mat_double::resize(int ir, int ic)
{
  m_nrows = ir;
  m_ncols = ic;
  m_Data.resize(ir*ic, 0.0);
}

void mat_double::Einheitsmatrix(int ir, int ic)
{
  (*this).resize(ir, ic);
  for(int i = 0; i < ir; i++)
  {
    for(int j = 0; j < ic; j++)
    {
      if(i == j)
      {
        (*this)(i, j)= 1;
      }
    }
  }
}

mat_double mat_double::inverse(int ir, int ic)
{
  mat_double inverse;
  mat_double kopie_matrix =(*this);
  mat_double i;                            // i = Einheitsmatrix
  i.Einheitsmatrix(ir, ic);
  gaussj(kopie_matrix, i);
  inverse = kopie_matrix;
  return inverse;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
vec_double operator*(mat_double a, vec_double b)             // Zeile aus der Matrix * Vektor        // freie Funktion
{
  vec_double Vret;
	Try
    vec_double VZeile;                 // Zeile aus der Matrix
    if(a.nrows()== b.isize())
    {
      for(int ir = 0; ir < a.nrows(); ir++)
      {
        a.get_row(ir, VZeile, 0);
        double d = VZeile *b;
        Vret.push_back(d);
      }
    }
    else
    {
      nrerror("Fehler:vec_double operator*");
    }
  Catch
  return Vret;
}


mat_double operator*(mat_double a, mat_double b)
{
  mat_double ret;
  Try
  ret.multiplikation(a, b);
  //return ret;
  Catch
  return ret;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define Mat_IO_DP mat_double
#define Vec_INT vec_int
#define DP double
void SWAP(double& d1, double& d2)
{
  double d = d1;
  d1 = d2;
  d2 = d;
}

void gaussj(mat_double &a, mat_double &b)                                                       // "gauss jordan" Methode zur multiplikation zweier Matritzen
{
	Try
    if(a.nrows()== a.ncols()&& a.nrows()== b.nrows()&& b.ncols()> 0)
    {
      int i, icol, irow, j, k, l, ll;
      DP big, dum, pivinv;
      int n = a.nrows();
      int m = b.ncols();
      vec_int indxc(n), indxr(n), ipiv(n);
      irow = 0;
      icol = 0;
      for(j = 0; j < n; j++)ipiv[j] = 0;
      for(i = 0; i < n; i++){
        big = 0.0;
        for(j = 0; j < n; j++)
          if(ipiv[j] != 1)
        for(k = 0; k < n; k++){
          if(ipiv[k] == 0){
            if(fabs(a[j][k])>= big){
              big = fabs(a[j][k]);
              irow = j;
              icol = k;
            }
          }
        }
        ++(ipiv[icol]);
        if(irow != icol){
          for(l = 0; l < n; l++)SWAP(a[irow][l], a[icol][l]);
          for(l = 0; l < m; l++)SWAP(b[irow][l], b[icol][l]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if(a[icol][icol] == 0.0)nrerror("gaussj:Singular Matrix");
        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for(l = 0; l < n; l++)a[icol][l] *= pivinv;
        for(l = 0; l < m; l++)b[icol][l] *= pivinv;
        for(ll = 0; ll < n; ll++)
        if(ll != icol){
          dum = a[ll][icol];
          a[ll][icol] = 0.0;
          for(l = 0; l < n; l++)a[ll][l] -= a[icol][l] * dum;
          for(l = 0; l < m; l++)b[ll][l] -= b[icol][l] * dum;
        }
      }
      for(l = n - 1; l >= 0; l--){
        if(indxr[l] != indxc[l])
          for(k = 0; k < n; k++)
          SWAP(a[k][indxr[l]], a[k][indxc[l]]);
      }
      return;
    }
    else
    {
      nrerror("Fehler:gaussj Zeilen und Spalten passen nicht");
    }
  Catch
  nrerror("Fehler:gaussj");
}

vec_double mat_double::gaussj_mv(mat_double a, vec_double& b)                                   // Methode zur multiplikation einer Matrix mit einem Vector
{
	Try
    vec_double xx;
    mat_double m1(b, eSpalte);
    /*mat_double m1(b.isize(), 1);
    m1.set_col(0, b, 0);*/
      gaussj(a, m1);
    m1.get_col(0, xx, 0);
    return xx;
  Catch
  nrerror("Fehler:mat_double::gaussj_mv");
  return vec_double();
}

void mat_double::add_row(vec_double vec)                                                        // Zeile am Ende eine Matrix hinzufügen und mit einem Vektor füllen
{
	Try
    if(vec.isize()== m_nrows)
    {
      int iz = nrows();   //zeile
      int isp = ncols();   //spalte
      mat_double kopie_data(iz+1, isp);
      for(int i = 0; i < m_nrows; i++)
      {
        for(int j = 0; j < m_ncols; j++)
        {
          kopie_data(i, j)=(*this)(i, j);
        }
      }
      kopie_data.set_row(iz, vec, 0);
      (*this)= kopie_data;
      return;
    }
    else
    {
      nrerror("mat_double::add_row , Dimensionen passen nicht");
    }
  Catch
nrerror("Catch mat_double::add_row");
}
void mat_double::insert_row(vec_double vec, int eZeile)           // einen Vektor, an einer Stelle-x(ZEILE), in die Matrix einfügen und die Matrix somit erweitern
{
	Try
    if(vec.isize()== m_nrows)
    {
      int iz = nrows();   //zeile
      int isp = ncols();   //spalte
      mat_double kopie_data(iz+1, isp);
      for(int i = 0; i < m_nrows -(m_nrows + 1 - eZeile); i++)
      {
        for(int j = 0; j < m_ncols ; j++)
        {
          kopie_data(i, j)=(*this)(i, j);
        }
      }
      kopie_data.set_row(eZeile - 1, vec, 0);
      for(int i = eZeile - 1; i < m_nrows; i++)
      {
        for(int j = 0; j < m_ncols; j++)
        {
          kopie_data(i+1, j)=(*this)(i, j);
        }
      }
      (*this)= kopie_data;
      return;
    }
    else
    {
      nrerror("mat_double::insert_row, Dimensionen passen nicht");
    }
  Catch
nrerror("Catch mat_double::insert_row");
//	Try                                                                                          // in m_Data einen vektor als Zeile einfügen
//    if(vec.isize()== m_ncols)
//    {
//      vec_double kopie_data;
//      int iz = nrows();   //zeile
//      int isp = ncols();   //spalte
//      kopie_data.resize((iz + 1)*isp, 0.0);
//      int k = 0;
//      for(int i = 0; i < m_Data.isize(); i++)
//      {
//        if(k == m_ncols*eZeile)
//        {
//          k = k + isp;
//        }
//        kopie_data[k] = m_Data[i];
//        k++;
//      }
//      for(int i = 0; i < vec.isize(); i++)
//      {
//        kopie_data[m_ncols*eZeile + i] = vec[i];
//      }
//      m_Data.resize((iz + 1)*isp, 0.0);
//      m_Data = kopie_data;
//      m_nrows++;
//      return;
//    }
//    else
//    {
//      nrerror("mat_double::insert_row");
//    }
//  Catch
//nrerror("Catch mat_double::insert_row");
}
void mat_double::add_col(vec_double vec)                                                       // Spalte am Ende eine Matrix hinzufügen und mit einem Vektor füllen
{
	Try
    if(vec.isize()== m_ncols)
    {
      int iz = nrows();   //zeile
      int isp = ncols();   //spalte
      mat_double kopie_data(iz, isp + 1);
      for(int i = 0; i < m_nrows; i++)
      {
        for(int j = 0; j < m_ncols; j++)
        {
          kopie_data(i, j)=(*this)(i, j);
        }
      }
      kopie_data.set_col(isp,  vec, 0);
      (*this)= kopie_data;
      return;
    }
    else
    {
      nrerror("mat_double::add_col , Dimensionen passen nicht");
    }
  Catch
nrerror("Catch mat_double::add_col");
}
void mat_double::insert_col(vec_double vec, int eSpalte)         // einen Vektor, an einer Stelle-x(SPALTE), in die Matrix einfügen und die Matrix somit erweitern
{
	Try
    if(vec.isize()== m_ncols)
    {
      int iz = nrows();   //zeile
      int isp = ncols();   //spalte
      mat_double kopie_data(iz, isp + 1);
      for(int i = 0; i < m_nrows ; i++)
      {
        for(int j = 0; j < m_ncols -(m_ncols+1 - eSpalte); j++)
        {
          kopie_data(i, j)=(*this)(i, j);
        }
      }
      kopie_data.set_col(eSpalte-1, vec, 0);
      for(int i = 0; i < m_nrows; i++)
      {
        for(int j = eSpalte-1; j < m_ncols; j++)
        {
          kopie_data(i, j+1)=(*this)(i, j);
        }
      }
      (*this)= kopie_data;
      return;
    }
    else
    {
      nrerror("mat_double::insert_col, Dimensionen passen nicht");
    }
  Catch
nrerror("Catch mat_double::insert_col");
}


void mat_double::multiplikation(mat_double a, mat_double b)
{
  try{
    mat_double c(a.m_nrows, b.m_ncols);
    if (a.m_ncols != b.m_nrows)
    {
      nrerror("Dimensionen von Spalte und Zeile sind falsch");
    }
    else
    {
      for (int iz = 0; iz < c.m_nrows; iz++)
      {
        for (int is = 0; is < c.m_ncols; is++)
        {

          c(iz, is) = 0;     // Summe an stell iz, is berechnet

          for (int i = 0; i < a.m_ncols; i++)
          {
            c(iz, is) += a(iz, i) *  b(i, is);   // i < a.m_ncols   & i< b.m_nrows
          }

        }
      }
    }

    (*this) = c;
  }  catch (...){}
}







					//////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////

nr_Optimierung::nr_Optimierung(void)
{
  m_anzahlIterationen = 50000;
  phi = (1.0 + sqrt(5.0)) / 2.0;
  resphi = 2 - phi;
}



          bool nr_Optimierung::GoldenerSchnitt_Start(nr_Functor1& f, double& x1, double& x3, int N)
          {
            double t, w, tmin = 0.0;
            double wmin = f(x1)+ 1.0;
            double ds =(x3 - x1)/ N;
            for(int i = 0; i <= N; i++)
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
              x1 = tmin - dds;
              x3 = tmin + dds;
              if(f(x1)> wmin && f(x3)> wmin)
              {
                return true;
              }
              dds += 0.1*ds;
            }
            return false;
          }
					//
          double nr_Optimierung::GoldenerSchnitt(nr_Functor1& f, double& x1, double& x3, double eps)
          {
            double w2, w4, x5, w5;
            if(fabs(x1 - x3)< eps)
            {
              return(x1 + x3)/ 2.0;
            }
            double x2 = x1 + resphi *(x3 - x1);
            w2 = f(x2);
						// Create a new possible center in the area between x2 and x3, closer to x2
            double x4 = x2 + resphi *(x3 - x2);
            w4 = f(x4);
            while(fabs(x1 - x3)> eps)
            {
              if(w4 < w2)
              {
                x1 = x2;
                x2 = x4;
                w2 = w4;
                if((x2 - x1)<(x3 - x2))
                {
                  x4 = x2 + resphi *(x3 - x2);
                  w4 = f(x4);
                }
                else
                {
                  x5 = x2 - resphi *(x2 - x1);
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
                if((x2 - x1)<(x3 - x2))
                {
                  x4 = x2 + resphi *(x3 - x2);
                  w4 = f(x4);
                }
                else
                {
                  x5 = x2 - resphi *(x2 - x1);
                  w5 = f(x5);
                  x4 = x2;
                  w4 = w2;
                  x2 = x5;
                  w2 = w5;
                }
              }
            }
            return(x1 + x3)/ 2.0;
          }
					//
					//
          double nr_Optimierung::GoldenerSchnittR(nr_Functor1& f, double& x1, double& x2, double& x3, double eps)
          {
            if(fabs(x1 - x3)< eps)
            {
              return(x1 + x3)/ 2.0;
            }
						// Create a new possible center in the area between x2 and x3, closer to x2
            double x4 = x2 + resphi *(x3 - x2);
            if(f.w(x4)< f.w(x2))
            {
              return GoldenerSchnittR(f, x2, x4, x3, eps);
            }
            else
            {
              return GoldenerSchnittR(f, x4, x2, x1, eps);
            }
            return(x1 + x3)/ 2.0;
          }


          void nr_Optimierung::make_amoeba_StartMatrix(mat_double & m,
            vec_double &p0,
            int NDim, double lambda)
          {
            m.set_rc(NDim + 1, NDim);
            m.set_row(0, p0, NDim);
            vec_double p;
            p.resize(p0.size());
            for (int i = 1; i <= NDim; i++)
            {
              p = p0;
              p.add(i - 1, lambda);
              m.set_row(i, p, NDim);
            }
          }


          void nr_Optimierung::get_psum(mat_double &p, vec_double &psum)
          {
            int i, j;
            double sum;

            int mpts = p.nrows();
            int ndim = p.ncols();
            for (j = 0; j<ndim; j++)
            {
              for (sum = 0.0, i = 0; i<mpts; i++)
              {
                sum += p(i, j);
              }
              psum.setat(j, sum);
            }
          }

          bool nr_Optimierung::nr_amoeba(int NDim, double ftol, double lambda,
            vec_double start,
            nr_FunctorN& f,
            vec_double& erg)
          {
            if (nr_amoeba_xmal(NDim, ftol, lambda, start, f, erg))
            {
              start = erg;
              if (nr_amoeba_xmal(NDim, ftol, lambda, start, f, erg))
              {
                lambda = 0.1*lambda;
                start = erg;
                if (nr_amoeba_xmal(NDim, ftol, lambda, start, f, erg))
                {
                  return true;
                }
              }
            }
            return false;
          }

          bool nr_Optimierung::nr_amoeba_xmal(int NDim, double ftol, double lambda,
            vec_double start,
            nr_FunctorN& f,
            vec_double& erg)
          {
            mat_double m;
            make_amoeba_StartMatrix(m, start, NDim, lambda);
            vec_double y(NDim + 1);
            vec_double pt(NDim);
            for (int i = 0; i <= NDim; i++)
            {
              m.get_row(i, pt, NDim);
              double dff = f(pt);
              y.setat(i, dff);
            }

            int nfunk;
            nr_amoeba(m, y, ftol, f, nfunk);
            double test;
            if (nfunk < m_anzahlIterationen)
            {
              erg = vec_double(NDim);
              for (int i = 0; i <= NDim; i++)
              {
                m.get_row(i, pt, NDim);
                test = f(pt);
                erg.add_vec(NDim, pt);
              }
              for (int i = 0; i<NDim; i++)
              {
                double d = erg.getat(i);
                d = d / (NDim + 1);
                erg.setat(i, d);
              }
              test = f(erg);
              return true;
            }
            return false;
          }

          void nr_Optimierung::nr_amoeba(mat_double &p, vec_double &y,
            const double ftol, nr_FunctorN& funk,
            int &nfunk)
          {
            const int NMAX = 5000;
            const double TINY = 1.0e-10;
            int i, ihi, ilo, inhi, j;
            double rtol, ysave, ytry;

            int mpts = p.nrows();
            int ndim = p.ncols();
            vec_double psum(ndim);
            nfunk = 0;
            get_psum(p, psum);
            for (;;) {
              ilo = 0;
              //ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);


              if (y[0] > y[1])
              {
                inhi = 1;
                ihi = 0;
              }
              else
              {
                inhi = 0;
                ihi = 1;
              }

              for (i = 0; i<mpts; i++) {
                if (y[i] <= y[ilo]) ilo = i;
                if (y[i] > y[ihi]) {
                  inhi = ihi;
                  ihi = i;
                }
                else if (y[i] > y[inhi] && i != ihi) inhi = i;
              }
              rtol = 2.0*fabs(y[ihi] - y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]) + TINY);
              if (rtol < ftol)
              {
                //SWAP(y[0],y[ilo]);
                y.swap(0, ilo);
                for (i = 0; i<ndim; i++) p.swap(0, i, ilo, i);//SWAP(p[0][i],p[ilo][i]);
                break;
              }
              if (nfunk >= NMAX)
              {
                nrerror("NMAX exceeded");
                return;
              }
              nfunk += 2;
              ytry = nr_amotry(p, y, psum, funk, ihi, -1.0);
              if (ytry <= y[ilo])
                ytry = nr_amotry(p, y, psum, funk, ihi, 2.0);
              else if (ytry >= y[inhi]) {
                ysave = y[ihi];
                ytry = nr_amotry(p, y, psum, funk, ihi, 0.5);
                if (ytry >= ysave) {
                  for (i = 0; i<mpts; i++) {
                    if (i != ilo) {
                      for (j = 0; j<ndim; j++)
                      {
                        double w = 0.5*(p(i, j) + p(ilo, j));
                        p.set(i, j, w);
                        psum.setat(j, w);
                      }
                      y[i] = funk(psum);
                    }
                  }
                  nfunk += ndim;
                  get_psum(p, psum);
                }
              }
              else --nfunk;
            }
          }



          double  nr_Optimierung::nr_amotry(mat_double &p,
            vec_double &y,
            vec_double &psum, nr_FunctorN& funk,
            const int ihi, const double fac)
          {
            int j;
            double fac1, fac2, ytry;

            int ndim = p.ncols();
            vec_double ptry(ndim);
            fac1 = (1.0 - fac) / ndim;
            fac2 = fac1 - fac;
            for (j = 0; j<ndim; j++)
              ptry.setat(j, psum[j]*fac1 - p(ihi, j)*fac2);
            ytry = funk(ptry);
            if (ytry < y[ihi])
            {
              y.setat(ihi, ytry);
              for (j = 0; j<ndim; j++) {
                psum.add(j, ptry[j] - p(ihi, j));
                p.set(ihi, j, ptry[j]);
              }
            }
            return ytry;
          }

					//////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////
          DP nr_Optimierung::rtsec(nr_Functor1& func, const DP x1, const DP x2, const DP eps)             // Methode zur Nullstellenbestimmung
          {
						Try
              const int MAXIT = 30;
              int j;
              DP fl, f, dx, xl, rts;
              fl = func(x1);
              f = func(x2);
              if(fabs(fl)< fabs(f)){
                rts = x1;
                xl = x2;
                SWAP(fl, f);
              }
              else {
                xl = x1;
                rts = x2;
              }
              for(j = 0; j<MAXIT; j++){
                dx =(xl - rts)*f /(f - fl);
                xl = rts;
                fl = f;
                rts += dx;
                f = func(rts);
                if(fabs(dx)< eps || f == 0.0)return rts;
              }
              nrerror("Maximum number of iterations exceeded in rtsec");
              return 0.0;
            Catch
            return 0.0;
          }
          DP nr_Optimierung::rtbis(nr_Functor1& func, const DP x1, const DP x2, const DP eps)             // Methode zur Nullstellenbestimmung
          {
						Try
              const int JMAX = 40;
              int j;
              DP dx, f, fmid, xmid, rtb;
              f = func(x1);
              fmid = func(x2);
              if(f*fmid >= 0.0)nrerror("Root must be bracketed for bisection in rtbis");
              rtb = f < 0.0 ?(dx = x2 - x1, x1):(dx = x1 - x2, x2);
              for(j = 0; j<JMAX; j++){
                fmid = func(xmid = rtb +(dx *= 0.5));
                if(fmid <= 0.0)rtb = xmid;
                if(fabs(dx)< eps || fmid == 0.0)return rtb;
              }
              nrerror("Too many bisections in rtbis");
              return 0.0;
            Catch
            return 0.0;
          }
          bool nr_Optimierung::zbrac(nr_Functor1& func, DP &x1, DP &x2)                                   // Methode zur Grenzbestimmung
          {
						Try
              const int NTRY = 50;
              const DP FACTOR = 1.6;
              int j;
              DP f1, f2;
              if(x1 == x2)nrerror("Bad initial range in zbrac");
              f1 = func(x1);
              f2 = func(x2);
              for(j = 0; j<NTRY; j++){
                if(f1*f2 < 0.0)return true;
                if(fabs(f1)< fabs(f2))
                  f1 = func(x1 += FACTOR*(x1 - x2));
                else
                  f2 = func(x2 += FACTOR*(x2 - x1));
              }
              return false;
            Catch
            return false;
          }
          void nr_Optimierung::zbrak(nr_Functor1& func, const DP x1, const DP x2, const int n, vec_double &xb1, vec_double &xb2, int &nroot) // Methode zur Grenzbestimmung
          {
						Try
              int i;
              DP x, fp, fc, dx;
              int nb = xb1.size();
              nroot = 0;
              dx =(x2 - x1)/ n;
              fp = func(x = x1);
              for(i = 0; i<n; i++){
                fc = func(x += dx);
                if(fc*fp <= 0.0){
                  xb1[nroot] = x - dx;
                  xb2[nroot++] = x;
                  if(nroot == nb)return;
                }
                fp = fc;
              }
            Catch
          }
          DP nr_Optimierung::bis_brac(nr_Functor1& func, DP x1, DP x2, const DP eps)   // sucht mit 2 verschiedenen Methoden neuen/passenden Grenzbereiche und führt  "rtbis" aus
					// wenn keine neue Grenzbereiche gefunden werden, wird eine Fehlermeldung ausgegeben.
          {
						Try
              DP fl, f;
              double ret = 0.0;
              bool bGültig = true;
              fl = func(x1);
              f = func(x2);
              if(fl*f > 0.0)  // prüft ob es ein Vorzeichenwechsel gibt
              {
                zbrac(func, x1, x2);
              }
              if(func(x1)*func(x2)>0.0)
              {
                int nroot = 0;
                vec_double xb1, xb2;
                int n = 5;
                zbrak(func, x1, x2, n, xb1, xb2, nroot);
                if(func(x1)*func(x2)> 0.0)
                {
                  bGültig = false;
                }
              }
              if(bGültig)
              {
                ret = rtbis(func, x1, x2, eps);
              }
              else
              {
                nrerror("kein passender Grenzbereich gefunden");
              }
              return ret;
            Catch
            nrerror("nr_Optimierung::bis_brac");
            return 0.0;
          }


   



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int nr_FunctorN_MinKugel::icounter = 0;
          //////////////////////////////////////////////////////////////////////////////////////////
          bool Test_vec_ok1(std::vector < cv3>& v)
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
          void nr_FunctorN_MinKugel::reset_vec_cv3(vec_double& vec)
          {
            v_cv3.resize(NN);
            v_cv3[0].set_polar(0.0, 0.0);
            v_cv3[1].set_polar(vec[0], 0.0);
            for (int i = 2; i < NN; i++)
            {
              v_cv3[i].set_polar(vec[2 * i - 3], vec[2 * i - 2]);
            }
            Test_vec_ok1(v_cv3);
          }
          //
          double nr_FunctorN_MinKugel::mindist()
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
          void nr_FunctorN_MinKugel::initVec(int iArt, vec_double& vec)
          {

            double pi = 3.1415926535897932384626433832795;
            double pi2 = 2 * pi;

            switch (iArt)
            {
						case -1: //Oktagon
						{
							NN = 6;
							vec.resize(9,0.0);
							vec[0] = 0.0;

							vec[1] = 0.0;
							vec[2] = 0.5*pi;

							vec[3] = 0.0;
							vec[4] = 1.0*pi;

							vec[5] = 0.0;
							vec[6] = 1.5*pi;

							vec[7] = -0.25*pi;
							vec[8] = 0.0;

						}
						break;
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
              double hw_akt = pi / (iKreise + 1);
              int iPktaufKreis = (NN - 2 * iODD) / iKreise;
              int iMittelKreis = (iKreise + 1) / 2;
              int ipktGes = (iODD == 1 ? 2 : 1) + iPktaufKreis*iKreise;
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
          double nr_FunctorN_MinKugel::operator()(vec_double& vec)
          {
            reset_vec_cv3(vec);
            double sum = 0.0;
            double zzz = 0.0;

						switch (m_iOptiZiel012)
            {
            case 11:
              icounter--;
              if (icounter < 0)
              {
								m_iOptiZiel012 = 1;
              }
            case 0:
            default:
              for (int i = 0; i < NN; i++)
                for (int j = i + 1; j < NN; j++)
                {
                  double d = v_cv3[i].Dist2(v_cv3[j]);
                  sum += 1.0 / d;
                }
              break;
            case 1:
              for (int i = 0; i < NN; i++)
                for (int j = i + 1; j < NN; j++)
                {
                  double d = v_cv3[i].Dist(v_cv3[j]);
                  if (d < zzz)
                  {
                    zzz = d;
                  }
                }
              break;
						case 2:
						{
							double dmin = 123456789.0;
							for (int i = 0; i < NN; i++)
							{
								for (int j = i + 1; j < NN; j++)
								{
									double d = v_cv3[i].Dist2(v_cv3[j]);
									if (d < dmin){ dmin = d; }
									sum += 1.0 / d;
								}
							}
							sum = 0.01*sum + 1.0 / dmin;
						}
							break;
						}
            return sum-zzz;
          }
          //
          //





/////////////////////////////////////////////////////
          ////////////////////////////////////////////


