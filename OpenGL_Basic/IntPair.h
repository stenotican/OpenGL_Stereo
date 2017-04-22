#ifndef _INTPAIR_H_
#define _INTPAIR_H_

#include <map>
#include <vector>

class CIntPair
{
public:
  int i1,i2;
public:
  CIntPair(void);
public:
  CIntPair(int _i1, int _i2);
  virtual ~CIntPair(void);
public:
  bool operator<(const CIntPair& ip) const;
  bool operator==(const CIntPair& ip) const;
};

//-------------------------------------------
//-------------------------------------------
class CIntMatrix :public std::map<CIntPair, int>
{
public:

public:
void setat(int i1,int i2, int iValue);

bool find(int i1, int i2) const;
bool find(CIntPair) const;
bool find(CIntPair, int& iValue) const;

};
//
class CIntMessageMatrix: public CIntMatrix
{
public:
  static System::String^ MessageText(int iMessage);
};
//
//
//
//-------------------------------------------
class CDoubleVector :public std::map<int, double>
{
public:
 int m_dim; 
public:
CDoubleVector(){m_dim =0;}
CDoubleVector(int N);
CDoubleVector(const CDoubleVector& objectSrc);
void operator=(const CDoubleVector& objectSrc);

void setat(int i1, double dValue);

bool find(int i1) const;
bool find(int i1, double& dValue) const;
double operator()(int i) const;
double getat(int i) const;
void add(int i, double d_add);
void add(int NDim, const CDoubleVector& objectSrc);
void swap(int i, int j);
double sum(int NDim);
};//-------------------------------------------
class CDoubleMatrix :public std::map<CIntPair, double>
{
public:
	int m_nrows; 
	int nrows() {return m_nrows;}
	int m_ncols;
	int ncols() {return m_ncols;}
	void set_rc(int rows, int cols) {m_nrows = rows;m_ncols = cols;}
	void set_row(int ir, CDoubleVector vec, int Ndim); 
	void get_row(int ir, CDoubleVector& vec, int Ndim); 

public:
CDoubleMatrix(){m_nrows = 0;m_ncols=0;}
CDoubleMatrix(int rows, int cols);

double operator()(int i1,int i2) const;
void setat(int i1,int i2, double dValue);

bool find(int i1, int i2) const;
bool find(CIntPair) const;
bool find(CIntPair, double& iValue) const;
bool find(int i1, int i2, double& dValue) const;

double getat(int i1,int i2) const;
void swap(int i1, int i2, int j1, int j2);

};

///////////////////////////////////////////////////////////////////////////////////

#endif