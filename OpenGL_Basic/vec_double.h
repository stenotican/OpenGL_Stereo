#ifndef _VEC_DOUBLE_
#define _VEC_DOUBLE_

#include <vector>
#include "cv3.h"

//
class nr_Functor1;
class nr_FunctorN;
//


class vec_int : public std::vector < int >
{
public:

  vec_int();

  vec_int(int N);
  vec_int(const vec_int& objectSrc);
  vec_int& vec_int::operator= (const vec_int& objectSrc);

  int isize()const{ return int(size()); }
};

//
//
//
//
class vec_double : public std::vector<double>
{

public:
  vec_double();
  vec_double(const double&x);
  vec_double(const double&x, const double&y);
  vec_double(const double&x, const double&y, const double&z);
  vec_double(const double&x, const double&y, const double&z, const double&w);
  

  vec_double(int N);
  vec_double(const vec_double& objectSrc);
  vec_double& vec_double::operator= (const vec_double& objectSrc);

  int isize()const{ return int(size()); }
                                          
  double maximum()const;                                           
  double minimum()const;
  double skalarprodukt(const vec_double& a, const vec_double& b);   // zwei Vektoren::(x*x)+(y*y)+(z*z)...
  double norm()const;                                               // ein Vektor::Wurzel aus(x²+y²+z²...)
  void normieren(double länge = 1.0);                               // Vektor auf eine bestimmte Länge setzen
  double summe()const;
  double abstand(const vec_double&a)const;                          // Abstand zw. Vektoren bestimmen
  double winkel(const vec_double& a, const vec_double& b)const;     // Winkel zw. 2 Vektoren berechnen
  double getat(int i) const;
  void setat(int index, double iWert);
  void add(int i, double d_add);                                    // Addieren
  void add_vec(int NDim, const vec_double& v);
  void swap(int i, int j);                                          // Tausch
  
  
};

vec_double operator+(const vec_double& a, const vec_double& b);     // Addition
vec_double operator- (const vec_double& a, const vec_double& b);    // Subtraktion
double operator*(const vec_double& a, const vec_double& b);         // Skalarprodukt
vec_double operator*(const double a, const vec_double& b);          // ein double *vektor z.B 3.0*(1/2/3)
vec_double operator-(const vec_double& a);                          // Vorzeichenwechsel

//
//
//
//
class mat_double
{
public:

  enum Zeile_Spalte { eZeile, eSpalte };
  static double dummy;



protected:
  int m_nrows;            // Zeilen
  int m_ncols;            // Spalten
  vec_double m_Data;

public:

  int isize()const{ return int(m_Data.size()); }


  mat_double();                                                     // Leere Matrix anlegen, Anz. Zeilen und Spalten werden auf 0 initialisiert
  mat_double(int Zeilen, int Spalten);                              // Matrix anlegen mit Zeilen und Spalten
                                                                    // (Zeilen und Spalten der) Matrix mit 0.0 initialisieren

  mat_double(vec_double v, Zeile_Spalte eZS);

  double & operator()(int iZ, int iS);                              // ()Operator
  double* operator[](int iZ);                                       // [] Operator


  void init();                                                      // die Matrix mit Werten füllen


public:            // nu_Optimierung


  int nrows() { return m_nrows; }
  int ncols() { return m_ncols; }


  void set_rc(int rows, int cols) { m_nrows = rows; m_ncols = cols; m_Data.resize(rows*cols, 0.0); }

  void swap(int z1, int sp1, int z2, int sp2);
  
  void set(int ir, int ic, double wert);                            // ein Wert einfügen
  void set_row(int ir, vec_double vec, int Ndim);                   // den Vektor in einer Zeile ( der Matrix) durch einen Anderen ersetzen
  void set_col(int ic, vec_double vec, int Ndim);                   // den Vektor in einer Spalte ( der Matrix) durch einen Anderen ersetzen
  double& get(int ir, int ic);                                      // ein Wert auslesen
  void get_row(int ir, vec_double& vec, int Ndim);                  // einen Vektor aus der Matrix ( eine Zeile) auslesen können
  void get_col(int ic, vec_double& vec, int Ndim);                  // einen Vektor aus der Matrix ( eine Spalte) auslesen können
  void resize(int ir, int ic);
  void Einheitsmatrix(int ir, int ic);
  mat_double inverse(int ir, int ic);                               // A^-1 

  vec_double gaussj_mv(mat_double a, vec_double& b);                // "gauss jordan"  mit Matrix und Vektor benutzen

  void add_row(vec_double vec);                                     // Zeile am Ende eine Matrix hinzufügen und mit einem Vektor füllen
  void insert_row(vec_double vec, int eZeile);                      // einen Vektor, an einer Stelle-x (ZEILE), in die Matrix einfügen und die Matrix somit erweitern 
  void add_col(vec_double vec);                                     // Spalte am Ende eine Matrix hinzufügen und mit einem Vektor füllen
  void insert_col(vec_double vec, int eSpalte);                     // einen Vektor, an einer Stelle-x (SPALTE), in die Matrix einfügen und die Matrix somit erweitern

  
 void multiplikation(mat_double a, mat_double);

};

vec_double operator*(mat_double a,  vec_double b);                  // freie Funktion  // eine Zeile aus der Matrix mit einem Vektor multiplizieren
mat_double operator*(mat_double a, mat_double b);
void gaussj(mat_double& a, mat_double &b);                          // freie Funktion  // "gauss jordan" mit zwei Matrizen benutzen


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





#define DP double

class nr_Optimierung
{
public:
  double phi;
  double resphi;
  int m_anzahlIterationen;
 
public:
  nr_Optimierung(void);
  virtual ~nr_Optimierung(void) {}
public:
  static void nrerror(System::String^ s){}
public:

//Methoden, die aus der vorigen nu_Optimierung.h stammen

//                   bereits implementiert

  bool   GoldenerSchnitt_Start(nr_Functor1& f, double& x1, double& x3, int N = 10);         
  double GoldenerSchnittR(nr_Functor1& f, double& x1, double& x2, double& x3, double eps);  
  double GoldenerSchnitt(nr_Functor1& f, double& x1, double& x3, double eps);               
  //
  //
  //              
  //
  void make_amoeba_StartMatrix(mat_double & p,                                            
    vec_double &p0, int NDim, double lambda);
  void get_psum(mat_double &p, vec_double &psum);
  ////
  bool nr_amoeba(int NDim, double ftol, double lambda,                                     
    vec_double start,
    nr_FunctorN& f,
    vec_double& erg);

  ////
  bool nr_amoeba_xmal(int NDim, double ftol, double lambda,                                
    vec_double start,
    nr_FunctorN& f,
    vec_double& erg);

  ////
  void nr_amoeba(mat_double &p, vec_double &y,                                           
    const double ftol, nr_FunctorN& funk,
    int &nfunk);

  ////
  double  nr_amotry(mat_double &p,                                                       
    vec_double &y,
    vec_double &psum,
    nr_FunctorN& funk,
    const int ihi, const double fac);

//
//
//
  static DP rtbis(nr_Functor1& func, const DP x1, const DP x2, const DP eps);                // aus "numerical_recipes/cpp/recipes"(library) kopiert
  static bool zbrac(nr_Functor1& func, DP &x1, DP &x2);                                      // aus "numerical_recipes/cpp/recipes"(library) kopiert
  static void zbrak(nr_Functor1& func, const DP x1, const DP x2, const int n,                // aus "numerical_recipes/cpp/recipes"(library) kopiert
                    vec_double &xb1, vec_double &xb2, int &nroot);
  static DP rtsec(nr_Functor1& func, const DP x1, const DP x2, const DP eps);                // aus "numerical_recipes/cpp/recipes"(library) kopiert

  static DP nr_Optimierung::bis_brac(nr_Functor1& func,  DP x1,  DP x2, const DP eps);       // sucht mit 2 verschiedenen Methoden neuen/passenden Grenzbereiche und führt  "rtbis" aus 
                                                                                             // wenn keine neue Grenzbereiche gefunden werden, wird eine Fehlermeldung ausgegeben.
};

//
//
//
//
class nr_Functor1
{
public:


  nr_Functor1(void){}
  virtual ~nr_Functor1(void){}
public:
  virtual double w(double t){ t; return 0.0; }
  virtual double operator()(double t){ t; return 0.0; }

};

//
//
//
//
class nr_FunctorN
{
public:
  nr_FunctorN(void){}
  virtual ~nr_FunctorN(void){}
public:
  virtual double operator()(vec_double& vec){ vec; return 0.0; }
};
//
//
//
class nr_FunctorN_MinKugel :public nr_FunctorN
{
public:
  static int icounter;

public:
  int m_iOptiZiel012;
  int NN; //Anzahl kugeln
  std::vector < cv3> v_cv3;
public:
	nr_FunctorN_MinKugel(){ m_iOptiZiel012 = 1; }
  void initVec(int iArt, vec_double& vec);
  void reset_vec_cv3(vec_double& vec);
  virtual double operator()(vec_double& vec);
  double mindist();
};
//
//
//
//
class nr_Functor1_sqrt:public nr_Functor1
{
public:
  int m_iZähler;                                                                             // wie oft eine Funktion zur Berechnung ausgeführt/aufgerufen wird  
  int m_iArt;                                                                                // welche Funktionsgleichung (case ) genommen wird
  nr_Functor1_sqrt(){ m_iZähler = 0; }                                                
  virtual double operator()(double t);

};
//



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////Dialog aufbauen//////////////////////////////////////////////////////////////////



#endif