#include "StdAfx.h"

#include "glm\glm.hpp"

#include "glm\vec2.hpp"
#include "glm\vec3.hpp"
#include "glm\vec4.hpp"
#include "glm\mat4x4.hpp"
#include "glm\gtc\matrix_transform.hpp"

#include "OpenGL_Base.h"



#include "nu_optimierung.h"
#include "MyData.h"

#include <math.h>

#include "vec_double.h"


using namespace OpenGLForm;


v3::v3(float r, float w)
{
	float t2 = 1 + r*r;
	float t1 = sqrt(t2);
	x = t1*cos(w);
	y = t1*sin(w);
	z = r / t1;
}


v3::v3(int iArt, float r, float w) //xyz-Normale
{
	float t2 = 1 + r*r;
	  x = -cos(w)*(t2);
		y = -sin(w)*(t2);
		z = r;
		float n = Norm();
		x = x / n;
		y = y / n;
		z = z / n;
}
//
void v3::set_wurm(float r, float w) //xyz-Wurm
{
	float t2 = 1 + r*r;
	float t1 = sqrt(t2);
	x = t1*cos(w);
	y = t1*sin(w);
	z = r / t1;
}
//
void v3::set_polar(double gw, double hw, double r)
{
	z = sin(hw);
	double xy = cos(hw);
	x = r*cos(gw)*xy;
	y = r*sin(gw)*xy;
	z = r*z;
}

//
void v3::set_p(v3^ vp)
{
	x = vp->x;
	y = vp->y;
	z = vp->z;
}

void v3::set_r(v3% vr)
{
	x = vr.x;
	y = vr.y;
	z = vr.z;
}

//
float v3::Dist(v3% v)
{
	return sqrt((v.x - x)*(v.x - x) + (v.y - y)*(v.y - y) + (v.z - z)*(v.z - z));
}
//
float v3::Norm()
{ 
	return sqrt(x*x + y*y + z*z); 
}
void v3::Normieren(float laenge)
{
	float nn = Norm();
	if (nn > 0.0)
	{
		x = x*laenge / nn;
		y = y*laenge / nn;
		z = z*laenge / nn;
	}
}


v3^ v3::langs_this(const v3% a)
{
	v3^ v1 = gcnew v3(*this);
	v1->Normieren(1.0);
	v3% b = *v1;
	float ff = a.x*b.x + a.y*b.y + a.z*b.z;
	v1->Normieren(ff);
	return v1;
}

v3^ v3::quer_this(const v3% a)
{
	v3^ v1 = langs_this(a);
	v3^ v2 = gcnew v3(a);
  
	v3% v1r = *v1;
	v3% v2r = *v2;

	return(v2r - v1r);
}

v3^ v3::operator+(const v3% b)
{
	v3^v = gcnew v3();
	v->x = x + b.x;
	v->y = y + b.y;
	v->z = z + b.z;
	return v;
}

v3^ v3::operator-(const v3% b)
{
	v3^v = gcnew v3();
	v->x = x - b.x;
	v->y = y - b.y;
	v->z = z - b.z;
	return v;
}

System::String^ v3::ToString()
{
	using namespace System::Globalization;
	CultureInfo^ MyCI = gcnew CultureInfo("en-US", false);
	NumberFormatInfo^ nfi = MyCI->NumberFormat;
	nfi->NumberDecimalSeparator = ".";
	return "<" + x.ToString(nfi) + " , " + y.ToString(nfi) + " , " + z.ToString(nfi) + ">";
}
//----------------------------------------------------------------
void listv3::add(String^ s3)
{
	array<String^>^ a = s3->Split();
	int is = a->GetLength(0);
	double x = is>0 ? System::Convert::ToDouble(a[0]) : 0.0;
	double y = is>0 ? System::Convert::ToDouble(a[1]) : 0.0;
	double z = is>0 ? System::Convert::ToDouble(a[2]) : 0.0;
	v3^ v = gcnew v3(x, y, z);
	Add(v);
}
//----------------------------------------------------------------
void listi3::add(String^ s3)
{
	array<String^>^ a = s3->Split();
	int is = a->GetLength(0);
	int x = is>0 ? System::Convert::ToInt32(a[0]) : 0.0;
	int y = is>0 ? System::Convert::ToInt32(a[1]) : 0.0;
	int z = is>0 ? System::Convert::ToInt32(a[2]) : 0.0;
	i3^ v = gcnew i3(x, y, z);
	Add(v);
}


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

v3^ OpenGLForm::operator*(const float a, const v3% b)
{
	v3^v = gcnew v3();
	v->x = a * b.x;
	v->y = a * b.y;
	v->z = a * b.z;
	return v;
}


float OpenGLForm::operator*(const v3% a, const v3% b)
{
	float f = a.x*b.x + a.y*b.y + a.z*b.z;
	return f;
}

v3^ OpenGLForm::operator%(const v3% a, const v3% b)
{
	v3^v = gcnew v3();
	v->x = a.y*b.z - a.z*b.y;
	v->z = a.z*b.x - a.x*b.z;
	v->x = a.x*b.y - a.y*b.x;
	return v;
}
//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
/*
float operator*(const v3^ a, const v3^ b)
{
	float f = a->x*b->x + a->y*b->y + a->z*b->z;
	return f;
}
*/
//
void Wurmloch::get_rw(v3^ s2, float% r, float% w, float% dist, v3^% s3)
{
	float x = s2->x;
	float y = s2->y;
	float z = s2->z;

	w = atan2(y, x);
	r = x*x + y*y - 1.0;

	if (s2->z > 0)
		r = sqrt(r);
	else
		r = -sqrt(r);

	s3 = gcnew v3(r, w);
	dist = s3->Dist(*s2);

}

double Wurmloch::PolynomWert(double x, std::vector<double>& ai)
{
  double w = 0.0;
  for (int i = ai.size() - 1; i >= 0; i--)
  {
    w = x*w + ai[i];
  }
  return w;
}

void Wurmloch::get_rw_polynom(double t01, double r1, double w1, double r2, double w2, std::vector<double>& ai, double& w, double& r)
{
  w = w1 + t01*(w2 - w1);
  double p = PolynomWert(t01, ai);

}


List<v3^>^ Wurmloch::GetLine1(int N, float r1, float w1, float r2, float w2)
	{
		List<v3^>^arr = gcnew List < v3^ >;
		for (int i = 0; i <= N; i++)
		{
			float t = float(i) / float(N);
			float r = r1 + t*(r2 - r1);
			float w = w1 + t*(w2 - w1);
			arr->Add(gcnew v3(r, w));

		}
		return arr;
	}
//
List<v3^>^ Wurmloch::GetLine2(int N, float r1, float w1, float r2, float w2)
{
	List<v3^>^arr = gcnew List < v3^ >;

	for (int i = 0; i <= N; i++)
	{
		float t = float(i) / float(N);
		float r = r1 + t*(r2 - r1);
		float w = w1 + t*(w2 - w1);
		arr->Add(gcnew v3(r, w));

	}
	return arr;
}

List<v3^>^ Wurmloch::GetLine3(int N, float r1, float w1, float r2, float w2)
{
	std::vector<double> arr_tw;
	GetArr_TW(N, r1, w1, r2, w2, arr_tw);

	List<v3^>^arr = MakeLine(arr_tw);
	return arr;
}
List<v3^>^ Wurmloch::MakeLine(std::vector<double>& arr_tw)
{
	List<v3^>^arr = gcnew List < v3^ >;
	for (int i = 0; i < arr_tw.size(); i+=2)
	{
		arr->Add(gcnew v3(arr_tw[i], arr_tw[i + 1]));
	}
	return arr;
}


void Wurmloch::GetArr_TW(int N, float r1, float w1, float r2, float w2, std::vector<double>& arr_tw)
{
	CDoubleVector vec_start(2 * N);

	for (int i = 1; i < 2*N; i+=2)
	{
		float t = float(i) / float(N);
		float r = r1 + t*(r2 - r1);
		float w = w1 + t*(w2 - w1);
		vec_start.setat(i-1, r);
		vec_start.setat(i, w);
	}

	nu_FunctorN_WurmlochTW fn;
	fn.t1 = r1;
	fn.t2 = r2;
	fn.w1 = w1;
	fn.w2 = w2;
	fn.iNZwischen = vec_start.size() / 2;

	nu_Optimierung opt;
	int NDim = vec_start.size();
	double ftol = 0.0001;
	double lambda = 0.01;
	CDoubleVector  vec_erg(2 * N);

	opt.nu_amoeba(NDim, ftol, lambda, vec_start, fn, vec_erg);

	arr_tw.clear();
	arr_tw.push_back(r1);
	arr_tw.push_back(w1);
	for (int i = 0; i < vec_start.size(); i++)
	{
		arr_tw.push_back(vec_erg[i]);
	}
	arr_tw.push_back(r2);
	arr_tw.push_back(w2);




}
//
//
List<v3^>^ Wurmloch::GetLine4Winkel(int N, float r1, float w1, float r2, float w2)
{
	CDoubleVector vec_start(N);

	nu_FunctorN_Wurmloch fn;
	fn.t1 = r1;
	fn.t2 = r2;
	fn.w1 = w1;
	fn.w2 = w2;
	fn.N = N;
	fn.initVec(0,vec_start);

	nu_Optimierung opt;
	int NDim = vec_start.size();
	double ftol = 0.0001;
	double lambda = 0.01;
	CDoubleVector  vec_erg(N);

	opt.nu_amoeba(NDim, ftol, lambda, vec_start, fn, vec_erg);

	List<v3^>^arr = gcnew List < v3^ >;

	for (int i = -1; i <= N; i++)
	{
		float r = fn.getT(i);
		float w = fn.getW(i, vec_erg);
		arr->Add(gcnew v3(r, w));

	}
	return arr;

}
//
List<v3^>^ Wurmloch::GetLine5_steps(int iArt, int N, float r1, float w1, float r2, float w2)
{
	if (iArt==2)
	{
		w2 = w2 - 6.283185307179586476925286766559; //=2 *3.1415926535897932384626433832795
	}


	std::vector<double> va(N+1);
	std::vector<double> vt(N+1);

	std::vector<cv3> vp(N + 1);

	for (int i = 0; i <=N; i ++)
	{
		double t = float(i) / float(N);
		double r = r1 + t*(r2 - r1);
		double w = w1 + t*(w2 - w1);
		va[i] = w;
		vt[i] = r;
		cv3 p(r, w);
		vp[i] = p;
	}


	int icount = 9000;
	while (icount > 0)
	{
		icount--;
		for (int i = 1; i < N; i++)
		{
			cv3 p1 = vp[i - 1];
			cv3 p2 = vp[i + 1];
			cv3 pi = vp[i];

			double d1a = p1.Dist(pi);
			double d2a = p2.Dist(pi);

			double db = max(d1a, d2a);

			cv3 m = 0.5*(p1+p2);
			cv3 delta = m - pi;
			cv3 df_dt(1, vt[i], va[i]);
			cv3 df_da(2, vt[i], va[i]);
			double delta_t = delta*df_dt;
			double delta_a = delta*df_da;

			double v_t = df_dt*df_dt;
			double v_a = df_da*df_da;
			if (v_t > 1.0E-10)
				delta_t = delta_t / v_t;
			else
				delta_t = 0.0;
			if (v_a > 1.0E-10)
				delta_a = delta_a / v_a;
			else
				delta_a = 0.0;

			cv3 pi_neu(vt[i] + delta_t, va[i] + delta_a);

			double d1b = p1.Dist(pi_neu);
			double d2b = p2.Dist(pi_neu);
			double dd2 = max(d1b, d2b);


			int ic = 20;
			while (ic > 0 && dd2 > db)
			{
				ic--;
				delta_t = 0.5*delta_t;
				delta_a = 0.5*delta_a;
				pi_neu.set(vt[i] + delta_t, va[i] + delta_a);
				d1b = p1.Dist(pi_neu);
				d2b = p2.Dist(pi_neu);
				dd2 = max(d1b, d2b);
			}

			vt[i] = vt[i] + delta_t;
			va[i] = va[i] + delta_a;
			pi.set(vt[i], va[i]);
			vp[i] = pi;
		}
	}


	List<v3^>^arr = gcnew List < v3^ >;

	for (int i = 0; i <= N; i++)
	{
		arr->Add(gcnew v3(vt[i], va[i]));
	}
	return arr;
}
//

//
List<v3^>^ Wurmloch::MinKugeln(int iArtStart, int Ni, int iAnzIter, int iOptiZiel) //kugeln auf einheitskugel
{
  

	nr_FunctorN_MinKugel fn;
  if (iOptiZiel == 1)
  {
    iOptiZiel = 2;
    fn.icounter = iAnzIter / 10;
  }
	fn.m_iOptiZiel012 = iOptiZiel;
	fn.NN = Ni;
	vec_double vec_start(2 * fn.NN - 3);
	fn.initVec(iArtStart,vec_start);
		nr_Optimierung opt;
		opt.m_anzahlIterationen = iAnzIter;
		int NDim = vec_start.size();
		double ftol = 1.0E-10;
		double lambda = 0.1;
		vec_double  vec_erg(NDim);
		vec_erg = vec_start;

		int iokk = 0;
		if (iAnzIter > 0)
		{
			if (opt.nr_amoeba(NDim, ftol, lambda, vec_start, fn, vec_erg))
			{
				iokk = 1;
			}
		}
	List<v3^>^arr = gcnew List < v3^ >;

	fn.reset_vec_cv3(vec_erg);
	MinDist = fn.mindist();

	for (int i = 0; i < fn.NN; i++)
	{
		arr->Add(gcnew v3(fn.v_cv3[i].x, fn.v_cv3[i].y, fn.v_cv3[i].z));
	}

	// Write the string array to a new file named "Kugeln(N).txt".
	String^ filename;
	filename = "KUGELN_" + fn.NN.ToString() + ".txt";
	IO::StreamWriter^ outputFile = gcnew IO::StreamWriter(filename);

	for (int i = 0; i < fn.NN; i++)
	{
		System::Globalization::CultureInfo^ ui = System::Globalization::CultureInfo::InvariantCulture;
		String^ line = "v " + fn.v_cv3[i].x.ToString("F4", ui) + " " + fn.v_cv3[i].y.ToString("F4", ui) + " " + fn.v_cv3[i].z.ToString("F4", ui);
		//line->Format("v {0,12:F4} {1,12:F4} {2,12:F4]", fn.v_cv3[i].x, fn.v_cv3[i].y, fn.v_cv3[i].z);
		outputFile->WriteLine(line);
	}

	//faces
	int N = fn.NN;
	outputFile->WriteLine(" ");
	for (int i1 = 0; i1 < N; i1++)
	{
		for (int i2 = i1 + 1; i2 < N; i2++)
		{
			for (int i3 = i2 + 1; i3 < N;i3++)
			{ 
				cebene e(fn.v_cv3[i1], fn.v_cv3[i2], fn.v_cv3[i3], true);
				bool ebene_ok = true;
				double dmin = 1.0E10;
				double dmax = -1.0E10;

				for (int i4 = 0; ebene_ok && (i4 < N); i4++)
				{
					if ((i4 != i1) && (i4 != i2) && (i4 != i3))
					{
						double dist = e.abstand(fn.v_cv3[i4]);
						if (dist > dmax)
						{
							dmax = dist;
						}
						if (dist < dmin)
						{
							dmin = dist;
						}
					}
				}

				if ((dmax > 0.001) && (dmin < -0.001))
				{
					ebene_ok = false;
				}


				if (ebene_ok)
				{
					String^ fline = "f " + (i1+1).ToString() + " " + (i2+1).ToString() + " " + (i3+1).ToString();
					outputFile->WriteLine(fline);
				}
			}
		}
	}


	outputFile->Close();

	return arr;

}

//
//
//
//
//int OpenGL_Base::iRender1or2 = 0;
//
OpenGL_Base::OpenGL_Base(System::Windows::Forms::Panel^  panel1, System::Windows::Forms::Panel^  panel2)
{
	iRender1or2 = 0;
	//
	//double a11=-2.0, a12=5.0, a21=3.0, a22=-4.0, b1=29.0, b2=-26.0, x1, x2;
	//solve2( a11,  a12,  a21,  a22,  b1,  b2,  x1,  x2);
  //
  m_iDoMouseMove = 0;
  m_iShowCubes_or_Gel=eCube;
	m_grav = new CGravitation();

  //
  rtri_angle=0.0;				// Angle for the triangle
	rquad_angle=0.0;				// Angle for the quad
	berg_angle=0.0;

  m_spinX = 0.0;
  m_spinY = 0.0;
  m_transX = 0.0;
	m_transYYY = 0.0;
	m_transZ = 0.0;

  //
	System::IntPtr pI1 = panel1->Handle;
	System::IntPtr pI2 = panel2->Handle;
	//m_hWnd = dynamic_cast<HWND>(panel->Handle.ToPointer());
	m_hWnd1 = (HWND)(panel1->Handle.ToPointer());
	m_hWnd2 = (HWND)(panel2->Handle.ToPointer());
	m_HDC1 = GetDC(m_hWnd1);
	m_HDC2 = GetDC(m_hWnd2);

	m_iSetPixelFormat1 = 0;
	m_iSetPixelFormat2 = 0;


	MyStartOpenGL1(panel1);
	MyStartOpenGL2(panel2);


}


OpenGL_Base::~OpenGL_Base(void)
{
}


void OpenGL_Base::MyStartOpenGL1(System::Windows::Forms::Panel^  panel)
{
	if (m_HDC1)
	{
		m_iSetPixelFormat1 = MySetPixelFormat1(m_HDC1);
		if (m_iSetPixelFormat1)
		{
			ReSizeGLScene(1,m_iSetPixelFormat1, panel->Width, panel->Height);
			InitGL();
		}
		if (m_iSetPixelFormat1 == 0)
		{
			m_iSetPixelFormat1 = MySetPixelFormat1(m_HDC1);
			if (m_iSetPixelFormat1)
			{
				ReSizeGLScene(1,m_iSetPixelFormat1, panel->Width, panel->Height);
				InitGL();
			}
		}
	}
}
void OpenGL_Base::MyStartOpenGL2(System::Windows::Forms::Panel^  panel)
{
	if (m_HDC2)
	{
		m_iSetPixelFormat2 = MySetPixelFormat2(m_HDC2);
		if (m_iSetPixelFormat2)
		{
			ReSizeGLScene(2,m_iSetPixelFormat2, panel->Width, panel->Height);
			InitGL();
		}
		if (m_iSetPixelFormat2 == 0)
		{
			m_iSetPixelFormat2 = MySetPixelFormat2(m_HDC2);
			if (m_iSetPixelFormat2)
			{
				ReSizeGLScene(2,m_iSetPixelFormat2, panel->Width, panel->Height);
				InitGL();
			}
		}
	}
}



GLint OpenGL_Base::MySetPixelFormat1(HDC hdc)
{
	static	PIXELFORMATDESCRIPTOR pfd =				// pfd Tells Windows How We Want Things To Be
	{
		sizeof(PIXELFORMATDESCRIPTOR),				// Size Of This Pixel Format Descriptor
		1,											// Version Number
		PFD_DRAW_TO_WINDOW |						// Format Must Support Window
		PFD_SUPPORT_OPENGL |						// Format Must Support OpenGL
		PFD_DOUBLEBUFFER,							// Must Support Double Buffering
		PFD_TYPE_RGBA,								// Request An RGBA Format
		16,										// Select Our Color Depth
		0, 0, 0, 0, 0, 0,							// Color Bits Ignored
		0,											// No Alpha Buffer
		0,											// Shift Bit Ignored
		0,											// No Accumulation Buffer
		0, 0, 0, 0,									// Accumulation Bits Ignored
		16,											// 16Bit Z-Buffer (Depth Buffer)  
		0,											// No Stencil Buffer
		0,											// No Auxiliary Buffer
		PFD_MAIN_PLANE,								// Main Drawing Layer
		0,											// Reserved
		0, 0, 0										// Layer Masks Ignored
	};

	GLint iPixelFormat = 0;

	// get the device context's best, available pixel format match 
	if ((iPixelFormat = ChoosePixelFormat(hdc, &pfd)) == 0)
	{
		MessageBox::Show("ChoosePixelFormat Failed");
		return 0;
	}

	// make that match the device context's current pixel format 
	if (SetPixelFormat(hdc, iPixelFormat, &pfd) == FALSE)
	{
		MessageBox::Show("SetPixelFormat Failed");
		return 0;
	}

	if (m_hglrc1 == NULL)
	{
		if ((m_hglrc1 = wglCreateContext(m_HDC1)) == NULL)
		{
			MessageBox::Show("wglCreateContext Failed");
			return 0;
		}
	}

	if ((wglMakeCurrent(m_HDC1, m_hglrc1)) == NULL)
	{
		MessageBox::Show("wglMakeCurrent Failed");
		return 0;
	}


	return 1;
}


GLint OpenGL_Base::MySetPixelFormat2(HDC hdc)
{
	static	PIXELFORMATDESCRIPTOR pfd =				// pfd Tells Windows How We Want Things To Be
	{
		sizeof(PIXELFORMATDESCRIPTOR),				// Size Of This Pixel Format Descriptor
		1,											// Version Number
		PFD_DRAW_TO_WINDOW |						// Format Must Support Window
		PFD_SUPPORT_OPENGL |						// Format Must Support OpenGL
		PFD_DOUBLEBUFFER,							// Must Support Double Buffering
		PFD_TYPE_RGBA,								// Request An RGBA Format
		16,										// Select Our Color Depth
		0, 0, 0, 0, 0, 0,							// Color Bits Ignored
		0,											// No Alpha Buffer
		0,											// Shift Bit Ignored
		0,											// No Accumulation Buffer
		0, 0, 0, 0,									// Accumulation Bits Ignored
		16,											// 16Bit Z-Buffer (Depth Buffer)  
		0,											// No Stencil Buffer
		0,											// No Auxiliary Buffer
		PFD_MAIN_PLANE,								// Main Drawing Layer
		0,											// Reserved
		0, 0, 0										// Layer Masks Ignored
	};

	GLint iPixelFormat = 0;

	// get the device context's best, available pixel format match 
	if ((iPixelFormat = ChoosePixelFormat(hdc, &pfd)) == 0)
	{
		MessageBox::Show("ChoosePixelFormat Failed");
		return 0;
	}

	// make that match the device context's current pixel format 
	if (SetPixelFormat(hdc, iPixelFormat, &pfd) == FALSE)
	{
		MessageBox::Show("SetPixelFormat Failed");
		return 0;
	}

	if (m_hglrc2 == NULL)
	{
		if ((m_hglrc2 = wglCreateContext(m_HDC2)) == NULL)
		{
			MessageBox::Show("wglCreateContext Failed");
			return 0;
		}
	}

	if ((wglMakeCurrent(m_HDC2, m_hglrc2)) == NULL)
	{
		MessageBox::Show("wglMakeCurrent Failed");
		return 0;
	}


	return 1;
}


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
		String^ OpenGL_Base::getString(int i)
		{
			v3^ v = (*arr_glob1)[i];
			return v->ToString();
		}

////////////////////////////////////////////////////////////////////////////////////////////
		glm::mat4 transform(
			glm::vec2 const & Orientation,
			glm::vec3 const & Translate,
			glm::vec3 const & Up)
		{
			glm::mat4 Projection = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.0f);
			glm::mat4 ViewTranslate = glm::translate(glm::mat4(1.0f), Translate);
			glm::mat4 ViewRotateX = glm::rotate(ViewTranslate, Orientation.y, Up);
			glm::mat4 View = glm::rotate(ViewRotateX, Orientation.x, Up);
			glm::mat4 Model = glm::mat4(1.0f);
			return Projection * View * Model;
		}
///////////////////////////////////////////////////
bool OpenGL_Base::InitGL(GLvoid)										// All setup for opengl goes here
		{
			glShadeModel(GL_SMOOTH);							// Enable smooth shading
			glClearColor(0.5f, 0.5f, 0.7f, 0.5f);				// Black background
			glClearDepth(1.0f);									// Depth buffer setup
			glEnable(GL_DEPTH_TEST);							// Enables depth testing
			glDepthFunc(GL_LEQUAL);								// The type of depth testing to do
			glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really nice perspective calculations



			{

			}

			return TRUE;										// Initialisation went ok


		}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

GLvoid OpenGL_Base::ReSizeGLScene(int iViewportNr, GLint iSetPixelformat, GLsizei width, GLsizei height)		// Resize and initialise the gl window
		{
			if (iSetPixelformat)
      {
			  if (height==0)										// Prevent A Divide By Zero By
			  {
				  height=1;										// Making Height Equal One
			  }

			  glViewport(0,0,width,height);						// Reset The Current Viewport

			  glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
			  glLoadIdentity();									// Reset The Projection Matrix

        

			  // Calculate The Aspect Ratio Of The Window
			  gluPerspective(45.0f,(GLfloat)width/(GLfloat)height,0.1f,100.0f);

				GLdouble delta_eye = 0.3333;

				switch (iViewportNr)
				{
				case 1:
					gluLookAt(-delta_eye, 0.0, 8.0, 0.0, 0.0, -0.0, 0.0, 1.0, 0.0);
					break;
				case 2:
					gluLookAt(delta_eye, 0.0, 8.0, 0.0, 0.0, -0.0, 0.0, 1.0, 0.0);
					break;
				}
        glRotated(m_spinX, 0.0, 1.0, 0.0);
        glRotated(m_spinY, 1.0, 0.0, 0.0);

			  glTranslated(m_transX, 0.0, 0.0);
				glTranslated(0.0, m_transYYY, 0.0);
				glTranslated(0.0, 0.0, m_transZ);

			  glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
			  glLoadIdentity();									// Reset The Modelview Matrix
        
		  }
    }
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
    float fkt(float angle, float x, float y)
    {
			int iArt = 1;
			switch(iArt)
			{
			case 1:
				{
					x = 3*x;
					y = 3*y;
					float z1 = 0.25*(2.0*cos(5.0*angle)/sqrt(1.0+0.3*(x*x+y*y)));
					float z2 = 0.25*(3.0*cos(1.4*angle+3.1415)/sqrt(1.0+0.3*((x+1)*(x+1)+(y-1)*(y-1))));
					float z3 = 0.25*(2.0*cos(2.333*angle+3.1415)/sqrt(1.0+0.3*((x-1)*(x-1)+(y+1)*(y+1))));
					float z4 = 0.25*(4.0*cos(7.666*angle+5.1415)/sqrt(1.0+0.3*((x-2)*(x-2)+(y+2)*(y+2))));
					return z1+z2+z3+z4;
					break;
				}
			case 2:
				{
					float d=1.5;
					float r1 = sqrt((x+d)*(x+d)+y*y);
					float z1 = 0.33*sin(-5.666*angle+3.0*r1)/(1.0+0.25*r1);
					float r2 = sqrt((x-d)*(x-d)+y*y);
					float z2 = 0.33*sin(-4.666*angle+3.0*r2)/(1.0+0.25*r2);
					return z1+z2;
				}
				break;
			}
			return 0.0;
    }

  void normale_fkt(float angle, float x, float y, float fkt(float , float , float), float h, float&nx,  float&ny,  float&nz)
  {
    float f1 = fkt(angle,x-h,y);
    float f2 = fkt(angle,x+h,y);
    float dx = 0.5*(f1+f2)/h;
     f1 = fkt(angle,x,y-h);
     f2 = fkt(angle,x,y+h);
    float dy = 0.5*(f1+f2)/h;

    float ll=sqrt(dx*dx+dy*dy+1.0);
    nx = dx/ll;
    ny=  -dy/ll;
    nz = 1.0/ll;

  }
		////////////////////////////////////////////////////////////////////////////////////////////
  void XglNormal3f(float x,float y, float z)
  {
    float t=sqrt(x*x+y*y+z+z);
    x=x/t;
    y=y/t;
    z=z/t;
    glNormal3f(x,y,z);

  }

  void X_normal_and_Vertex(float gr, float x,float y, float z)
  {
    float t=sqrt(x*x+y*y+z*z);
    glNormal3f(x/t,y/t,z/t);
    glVertex3f(gr*x,gr*y,gr*z);
  }

	void X_normal_and_Vertex(float gr, float x, float y, float z, cv3 p)
	{
		float t = sqrt(x*x + y*y + z*z);
		glNormal3f(x / t, y / t, z / t);
		glVertex3f(gr*x+p.x, gr*y+p.y, gr*z+p.z);
	}

	void X_Kugel(float gr, float radius, float x, float y, float z)
	{
			GLUquadricObj *quadric;
			quadric = gluNewQuadric();
			glPushMatrix();
			glTranslated(gr*x, gr*y, gr*z);
			glColor3f(1.0f, 1.0f, 1.0f);
			gluSphere(quadric, radius, 20, 20);

  			glPopMatrix();
		}

	void X_Vertex(float gr, float x, float y, float z)
	{
		glVertex3f(gr*x, gr*y, gr*z);
	}


	void X_Vertex(float gr, cv3 c)
	{
		glVertex3f(gr*c.x, gr*c.y, gr*c.z);
	}





	void Wurm_normal_and_Vertex(float gr, float r, float a)
	{
		float sq2 = (1.0f + r*r);
		float sq = sqrt(sq2);
		float x, y, z;
		x = sq*cos(a);
		y = sq*sin(a);
		z = r / sq;
		float nx, ny, nz;
		nx = -cos(a)/sq2;
		ny = -sin(a)/sq2;
		nz = r ;
		float n2 = nx*nx + ny*ny + nz*nz;
		n2 = sqrt(n2);
		nx = nx / n2;
		ny = ny / n2;
		nz = nz / n2;

		glNormal3f(nx,ny,nz);
		glVertex3f(x,y,z);

	}

	void Wurm_normal_Line(float gr, float r, float a)
	{
		float sq2 = (1.0f + r*r);
		float sq = sqrt(sq2);
		float x, y, z;
		x = sq*cos(a);
		y = sq*sin(a);
		z = r / sq;
		float nx, ny, nz;
		nx = -cos(a) / sq2;
		ny = -sin(a) / sq2;
		nz = r;
		float n2 = nx*nx + ny*ny + nz*nz;
		n2 = sqrt(n2);
		nx = nx / n2;
		ny = ny / n2;
		nz = nz / n2;

		nx = -0.2*nx + x;
		ny = -0.2*ny + y;
		nz = -0.2*nz + z;

		glVertex3f(x, y, z);
		glVertex3f(nx, ny, nz);

	}


	////////////////////////////////////////////////////////////////////////////////////////////
	System::Void OpenGL_Base::Render_Planeten(System::Void)
	{
		
		
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); //clear with background color
			glEnable(GL_LIGHTING); //enable the lighting system
			{ glEnable(GL_LIGHT0); //enable a light source. Other sources are GL_LIGHT1,..,GL_LIGHT7
			GLfloat DiffuseLight_color[] = { 1.0, 0.0, 0.0, 1.0 }; //Specify the color of light
			glLightfv(GL_LIGHT0, GL_DIFFUSE, DiffuseLight_color);
			GLfloat light_position[] = { -1.0, 1.0, -4.0, 1.0 }; //Specify the position of light source
			glLightfv(GL_LIGHT0, GL_POSITION, light_position);
			}
 glEnable(GL_LIGHT1); //enable 2nd light source
 GLfloat DiffuseLight_colorA[] = { 0.0, 0.0, 1.0, 1.0 }; //Specify the color of light
 GLfloat DiffuseLight_colorB[] = { 1.0, 0.0, 0.0, 1.0 }; //Specify the color of light
 GLfloat DiffuseLight_colorC[] = { 0.0, 1.0, 0.0, 1.0 }; //Specify the color of light
GLfloat light_position1[] = { 1.0, 1.0, -4.0, 1.0 }; //Specify 1st position of light source
GLfloat light_position2[] = { -1.0, 1.0, -4.0, 1.0 }; //Specify 2nd position of light source
if (true)
{
	glLightfv(GL_LIGHT1, GL_POSITION, light_position1);
}
else
{
	glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
}

glEnable(GL_COLOR_MATERIAL);
glColorMaterial(GL_FRONT, GL_DIFFUSE);
		
if (m_grav)
{

	glLightfv(GL_LIGHT1, GL_DIFFUSE, DiffuseLight_colorA);
	Render_Sphere(m_grav->p[0]);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, DiffuseLight_colorB);
	Render_Sphere(m_grav->p[1]);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, DiffuseLight_colorC);
	Render_Sphere(m_grav->p[2]);

	m_grav->Go();
}
	}
		////////////////////////////////////////////////////////////////////////////////////////////
	System::Void OpenGL_Base::Render_Sphere(cv3 p)
	{

		X_Kugel(1.0f, p.r, p.x, p.y, p.z);

		return;



		glBegin(GL_TRIANGLES);								// Start drawing a triangle


		float pi = 3.14159265359;
		float step1 = 0.25;
		float step2 = pi / 15;
		float groesse = 2.0;
		int ii = 0;
		for (float hw = -1.3; hw <= 1.3; hw += step1)
		{
			ii++;
			for (float gw = -0.0; gw <= 2.0*pi - step2 + 0.001; gw += step2)
			{
				ii++;
				if (ii % 2 == 0)
					glColor3f(1.0f, 0.0f, 0.0f);						// Red
				else
					glColor3f(0.0f, 0.0f, 1.0f);						// Blue

				float z1 = sin(hw);
				float xy1 = cos(hw);
				float x1 = xy1*cos(gw);
				float y1 = xy1*sin(gw);
				float x2 = xy1*cos(gw + step2);
				float y2 = xy1*sin(gw + step2);
				float z3 = sin(hw + step1);
				float xy3 = cos(hw + step1);
				float x3 = xy3*cos(gw);
				float y3 = xy3*sin(gw);
				float x4 = xy3*cos(gw + step2);
				float y4 = xy3*sin(gw + step2);

				X_normal_and_Vertex(groesse, x1, y1, z1);
				X_normal_and_Vertex(groesse, x2, y2, z1);
				X_normal_and_Vertex(groesse, x3, y3, z3);

				X_normal_and_Vertex(groesse, x2, y2, z1);
				X_normal_and_Vertex(groesse, x3, y3, z3);
				X_normal_and_Vertex(groesse, x4, y4, z3);
			}
		}
		glEnd();											// Done drawing the pyramid

	}
	////////////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////
		System::Void OpenGL_Base::Render_Sphere(System::Void)
	{

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
		glLoadIdentity();									// Reset the current modelview matrix

		glBegin(GL_TRIANGLES);								// Start drawing a triangle


		float pi = 3.14159265359;
		float step1 = 0.25;
		float step2 = pi / 15;
		float groesse = 2.0;
		int ii = 0;
		for (float hw = -1.3; hw <= 1.3; hw += step1)
		{
			ii++;
			for (float gw = -0.0; gw <= 2.0*pi - step2 + 0.001; gw += step2)
			{
				ii++;
				if (ii % 2 == 0)
					glColor3f(1.0f, 0.0f, 0.0f);						// Red
				else
					glColor3f(0.0f, 0.0f, 1.0f);						// Blue

				float z1 = sin(hw);
				float xy1 = cos(hw);
				float x1 = xy1*cos(gw);
				float y1 = xy1*sin(gw);
				float x2 = xy1*cos(gw + step2);
				float y2 = xy1*sin(gw + step2);
				float z3 = sin(hw + step1);
				float xy3 = cos(hw + step1);
				float x3 = xy3*cos(gw);
				float y3 = xy3*sin(gw);
				float x4 = xy3*cos(gw + step2);
				float y4 = xy3*sin(gw + step2);

				X_normal_and_Vertex(groesse, x1, y1, z1);
				X_normal_and_Vertex(groesse, x2, y2, z1);
				X_normal_and_Vertex(groesse, x3, y3, z3);

				X_normal_and_Vertex(groesse, x2, y2, z1);
				X_normal_and_Vertex(groesse, x3, y3, z3);
				X_normal_and_Vertex(groesse, x4, y4, z3);



			}
		}
		glEnd();											// Done drawing the pyramid

	}
	////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////
	System::Void OpenGL_Base::Render_Wurmloch(System::Void)
	{
		int iMat = 0;
		switch (iMat)
		{
		case 1:
		{
			// set up light
			float ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
			float diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
			float specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
			float position[] = { 200.0f, 300.0f, 100.0f, 1.0f };

			glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
			glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
			glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
			glLightfv(GL_LIGHT0, GL_POSITION, position);
			glEnable(GL_LIGHT0);
			glEnable(GL_LIGHTING);


			float mat_ambient1[] = { 0.129412f, 0.223529f, 0.027451f };
			float mat_diffuse1[] = { 0.780392f, 0.568627f, 0.113725f };
			float mat_specular1[] = { 0.992157f, 0.941176f, 0.807843f };
			float shine1 = 13.8f;

			//Emerald
			float mat_ambient2[] = { 0.0215f, 0.1745f, 0.9215f, 0.55f };
			float mat_ambient2a[] = { 0.0215f, 0.1745f, 0.0215f, 0.55f };
			float mat_diffuse2[] = { 0.07568f, 0.61424f, 0.07568f, 0.55f };
			float mat_specular2[] = { 0.633f, 0.727811f, 0.633f, 0.55f };
			float shine2 = 76.8f;

			glEnable(GL_COLOR_MATERIAL);
			glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);



			glShadeModel(GL_FLAT);


		}
		default:
			glDisable(GL_COLOR_MATERIAL);
		break;
		};



		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
		glLoadIdentity();									// Reset the current modelview matrix


		for (int iWasZeichnen = 1; iWasZeichnen <= 3; iWasZeichnen+=2)
		{
			switch (iWasZeichnen)
			{
			case 11:
			break;


			case 1: //fläche
			{
				glBegin(GL_TRIANGLES);								// Start drawing a triangle


				float pi = 3.14159265359;
				float step1 = 0.25;
				float step2 = pi / 15;
				float groesse = 2.0;
				int ii = 0;
				for (float rr = -5.3; rr <= 5.3; rr += step1)
				{
					ii++;
					for (float gw = -0.0; gw <= 2.0*pi - step2 + 0.001; gw += step2)
					{
						ii++;
						if (ii % 2 == 0)
							glColor3f(1.0f, 0.0f, 0.0f);						// Red
						else
							glColor3f(0.0f, 0.0f, 1.0f);						// Blue

						Wurm_normal_and_Vertex(groesse, rr, gw);
						Wurm_normal_and_Vertex(groesse, rr + step1, gw);
						Wurm_normal_and_Vertex(groesse, rr, gw + step2);

						Wurm_normal_and_Vertex(groesse, rr + step1, gw);
						Wurm_normal_and_Vertex(groesse, rr + step1, gw + step2);
						Wurm_normal_and_Vertex(groesse, rr, gw + step2);


					}
				}
				glEnd();											// Done drawing the pyramid
			}
			break;
			case 2: //normalen
			{
				glColor3f(1.0f, 1.0f, 0.0f);						// Yellow
				glLineWidth(2.0f);

				glBegin(GL_LINES);
				float pi = 3.14159265359;
				float step1 = 0.25;
				float step2 = pi / 15;
				float groesse = 2.0;
				int ii = 0;
				for (float rr = -5.3; rr <= 5.3; rr += step1)
				{
					ii++;
					for (float gw = -0.0; gw <= 2.0*pi - step2 + 0.001; gw += step2)
					{
						ii++;


						Wurm_normal_Line(groesse, rr, gw);
						Wurm_normal_Line(groesse, rr + step1, gw);
						Wurm_normal_Line(groesse, rr, gw + step2);

						Wurm_normal_Line(groesse, rr + step1, gw);
						Wurm_normal_Line(groesse, rr + step1, gw + step2);
						Wurm_normal_Line(groesse, rr, gw + step2);


					}
				}

				glEnd();
				break;
			}
			case 3:
			{
				int N = 75;
				if (arr_glob1 == nullptr)
				{
					arr_glob1 = Wurmloch::GetLine5_steps(1, N, 4.0f, 0.0f, 4.0f, 2.0f);
				}
				if (arr_glob2 == nullptr)
				{
					arr_glob2 = Wurmloch::GetLine5_steps(1,N, 4.0f, 0.0f, -4.0f, 2.0f);
				}
				if (arr_glob3 == nullptr)
				{
					arr_glob3 = Wurmloch::GetLine5_steps(2, N, 4.0f, 0.0f, -4.0f, 2.0f);
				}

				List<v3^>^ arr = arr_glob1;
				glColor3f(0.99f, 1.0f, 0.01f);						// Yellow
				glLineWidth(15.0f);
				glBegin(GL_LINES);
				for (int i = 0; i < N-1; i++)
				{
					v3^ a = arr[i];
					v3^ b = arr[i + 1];
					glVertex3f(a->x, a->y, a->z);
					glVertex3f(b->x, b->y, b->z);
				}
				glEnd();

				arr = arr_glob2;
				glColor3f(0.99f, 1.0f, 0.01f);						// Yellow
				glLineWidth(15.0f);
				glBegin(GL_LINES);
				for (int i = 0; i < N-1; i++)
				{
					v3^ a = arr[i];
					v3^ b = arr[i + 1];
					glVertex3f(a->x, a->y, a->z);
					glVertex3f(b->x, b->y, b->z);
				}
				glEnd();

				arr = arr_glob3;
				glColor3f(0.99f, 0.5f, 0.5f);						// Yellow
				glLineWidth(15.0f);
				glBegin(GL_LINES);
				for (int i = 0; i < N-1; i++)
				{
					v3^ a = arr[i];
					v3^ b = arr[i + 1];
					glVertex3f(a->x, a->y, a->z);
					glVertex3f(b->x, b->y, b->z);
				}
				glEnd();



			}
				break;

			case 4:
			{
				int N = 50;
				if (arr_glob1 == nullptr)
				{
					arr_glob1 = Wurmloch::GetLine4Winkel(N, 5.0f, 0.5f, -4.0f, 0.6);
				}
				List<v3^>^ arr = arr_glob1;
				glColor3f(0.99f, 0.88f, 0.05f);						// Yellow
				glLineWidth(15.0f);
				glBegin(GL_LINES);
				for (int i = 0; i < N; i++)
				{
					v3^ a = arr[i];
					v3^ b = arr[i + 1];
					glVertex3f(a->x, a->y, a->z);
					glVertex3f(b->x, b->y, b->z);
				}
			}
				glEnd();
				break;
			}
		}

	}
	//
	//
	System::Void OpenGL_Base::Render_MinKugeln(System::Void)
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
		glLoadIdentity();									// Reset the current modelview matrix
		{
			int N = m_data.m_iAnzahlKugeln;
			int iFormFaktor = m_data.m_iFormfaktor;
			if (arr_glob1 == nullptr)
			{
        arr_glob1 = Wurmloch::MinKugeln(iFormFaktor, N, m_data.m_iMaxIterations, m_data.m_iOptiZiel);
				m_data.m_dMinDist = Wurmloch::MinDist;
			}
			List<v3^>^ arr = arr_glob1;
			glColor3f(0.99f, 0.88f, 0.05f);						// Yellow
			glLineWidth(15.0f);
			glBegin(GL_LINES);
			for (int i = 0; i < N ; i++)
				for (int i2 = i + 1; i2 < N; i2++)
				{
					v3^ a = arr[i];
					v3^ b = arr[i2];
					double d = a->Dist(*b);
					if (d < 1.5*Wurmloch::MinDist)
					{
						glColor3f(0.09f, 0.98f, 0.05f);						// Green
						if (d < 1.25*Wurmloch::MinDist)
						{
							glColor3f(0.99f, 0.88f, 0.05f); //Yellow
						}


						if (d < 1.1*Wurmloch::MinDist)
						{
							glColor3f(0.99f, 0.08f, 0.05f); //red
						}

							X_Vertex(2.0, a->x, a->y, a->z);
							X_Vertex(2.0, b->x, b->y, b->z);
						
					}
				}
			glEnd();

			glEnable(GL_LIGHTING);
			glEnable(GL_LIGHT0);
			glEnable(GL_FLAT);
			for (int i = 0; i < N; i++)
			{
					v3^ a = arr[i];
					X_Kugel(2.0, 0.25, a->x, a->y, a->z);
			}
			glDisable(GL_LIGHTING);
		}


		cv3 p1(0.0,0,1);
		cv3 p2(0.0,0,-1);
		cv3 q1(10.0,0,0);
		cv3 q2(-10.0,0,0);
		double p, q;
		cv3 px, qx;
		verbindungskante(p1, p2, q1, q2, p, q, px, qx);
		glColor3f(0.99f, 0.88f, 0.95f);						// Yellow
		glLineWidth(10.0f);
		glBegin(GL_LINES);
		X_Vertex(2.0, p1);
		X_Vertex(2.0, p2);
		//X_Vertex(2.0, q1);
		//X_Vertex(2.0, q2);
		glColor3f(0.099f, 0.088f, 0.95f);						// Yellow
		//X_Vertex(2.0, px);
		//X_Vertex(2.0, qx);
		glEnd();


	}
	///////////////////////////////////////////////////////////////////////////////////////////
	void SetColor(cv3& nn)
	{
		cv3 n = nn;
		double dmax = nn.Max_xyz();
		if (dmax > 0.0)
		{
			n = (1.0 / dmax)*n;

		}
		glColor3d(0.5 + 0.5*n.x, 0.5 + 0.5*n.y, 0.5 + 0.5*n.z);
	}
	////////////////////////////////////////////////////////////////////////////////////////////
	/*
	http://people.sc.fsu.edu/~jburkardt/data/obj/obj.html
	*/
	System::Void OpenGL_Base::Render_Archimedes(System::Void)
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
		glLoadIdentity();									// Reset the current modelview matrix
		if (m_Obj)
		{
			double faktor = 2.0;
			int iStatus = 0;
			int i = 0;
			int iface = m_Obj->m_f[i];
			while (iface <0)
			{
				switch (iface)
				{
				case -3:
					glBegin(GL_TRIANGLES);								// Start drawing a triangle
					break;
				case -4:
					glBegin(GL_QUADS);								// Start drawing a triangle
					break;
				default:
					glBegin(GL_POLYGON);								// Start drawing a triangle
					break;
				}

				int j = i + 1;
				int k1 = m_Obj->m_f[j] - 1;
				int k2 = m_Obj->m_f[j + 1] - 1;
				int k3 = m_Obj->m_f[j + 2] - 1;
				cv3 nn = (m_Obj->m_v[k2] - m_Obj->m_v[k1]) % (m_Obj->m_v[k3] - m_Obj->m_v[k1]);
				nn.Normieren(1.0);
				//glColor3d(0.5 + 0.5*nn.x, 0.5 + 0.5*nn.y, 0.5 + 0.5*nn.z);
				SetColor(nn);

				for (j = i + 1; j <= i + abs(iface); j++)
				{
					/*
					switch (j % 4)
					{
					case 1:
						glColor3f(1.0f, 0.0f, 0.0f);
						break;
					case 2:
						glColor3f(1.0f, 1.0f, 0.0f);
						break;
					case 3:
						glColor3f(0.0f, 0.0f, 1.0f);
						break;
					case 0:
						glColor3f(0.0f, 1.0f, 0.0f);
						break;
					}
					*/
					int k = m_Obj->m_f[j];
					glVertex3d(faktor*m_Obj->m_v[k - 1].x, faktor*m_Obj->m_v[k - 1].y, faktor*m_Obj->m_v[k - 1].z);
				}
				glEnd();
				i = i + abs(iface) + 1;
				iface = 0;
				if (i < m_Obj->m_f.size())
				{
					iface = m_Obj->m_f[i];
				}
			}
		}
		else
		{
			glBegin(GL_TRIANGLES);								// Start drawing a triangle
			glColor3f(1.0f, 0.0f, 0.0f);						// Red
			glVertex3f(0.0f, 1.0f, 0.0f);					// Top Of triangle (front)
			glColor3f(0.0f, 1.0f, 0.0f);						// Green
			glVertex3f(-1.0f, -1.0f, 1.0f);					// Left of triangle (front)
			glColor3f(0.0f, 0.0f, 1.0f);						// Blue
			glVertex3f(1.0f, -1.0f, 1.0f);					// Right of triangle (front)
			glEnd();
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////
  		System::Void OpenGL_Base::Render_Sphere_Material(System::Void)
      {

   GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
   GLfloat mat_ambient[] = { 0.7, 0.7, 0.7, 1.0 };
   GLfloat mat_ambient_color[] = { 0.8, 0.8, 0.2, 1.0 };
   GLfloat mat_diffuse[] = { 0.1, 0.5, 0.8, 1.0 };
   GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
   GLfloat no_shininess[] = { 0.0 };
   GLfloat low_shininess[] = { 5.0 };
   GLfloat high_shininess[] = { 100.0 };
   GLfloat mat_emission[] = {0.3, 0.2, 0.2, 0.0};


      GLfloat light_ambient1[] = { 1.0, 1.0, 0.0, 1.0 };
      GLfloat light_diffuse[] = { 0.0, 0.0, 1.0, 1.0 };
      GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
      GLfloat light_position[] = { 2.0, 2.0, 5.0, 0.0 };
      GLfloat light_position1[] = { 5.0, -2.0, -5.0, 0.0 };

      GLfloat mat_shininess[] = { 50.0 };
      glClearColor (0.0, 0.0, 0.0, 0.0);
      glShadeModel (GL_SMOOTH);

      glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
      glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);


      //glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
      //glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
      glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);

      glLightfv(GL_LIGHT0, GL_POSITION, light_position);


      glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
      glLightfv(GL_LIGHT1, GL_POSITION, light_position1);

      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
      glEnable(GL_LIGHT1);
      glEnable(GL_DEPTH_TEST);

			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
			glLoadIdentity();									// Reset the current modelview matrix

			glBegin(GL_TRIANGLES);								// Start drawing a triangle


      float pi=3.14159265359;
			float step1 = 0.25;
			float step2 = pi/15;
			float groesse = 2.0;
			int ii = 0;
			for (float hw = -1.3; hw <= 1.3; hw += step1)
      {
        ii++;
				for (float gw = -0.0; gw <= 2.0*pi-step2+0.001; gw += step2)
        {
					ii++;
          if(ii%2)
          {
   glMaterialfv(GL_FRONT, GL_AMBIENT, no_mat);
   glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
   glMaterialfv(GL_FRONT, GL_SPECULAR, no_mat);
   glMaterialfv(GL_FRONT, GL_SHININESS, no_shininess);
   glMaterialfv(GL_FRONT, GL_EMISSION, no_mat); 
   glEnable(GL_COLOR_MATERIAL);
glColorMaterial(GL_FRONT, GL_DIFFUSE);
/* now glColor* changes diffuse reflection  */
glColor3f(0.2, 0.5, 0.8);
          }
          else
          {
   glMaterialfv(GL_FRONT, GL_AMBIENT, no_mat);
   glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
   glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);
   glMaterialfv(GL_FRONT, GL_EMISSION, no_mat);
   glEnable(GL_COLOR_MATERIAL);
glColorMaterial(GL_FRONT, GL_SPECULAR);
/* glColor* no longer changes diffuse reflection  */
/* now glColor* changes specular reflection  */
glColor3f(0.9, 0.0, 0.2);
          }

          float z1 = sin(hw);
          float xy1 = cos(hw);
          float x1 = xy1*cos(gw);
          float y1 = xy1*sin(gw);
          float x2 = xy1*cos(gw+step2);
          float y2 = xy1*sin(gw+step2);
          float z3 = sin(hw+step1);
          float xy3 = cos(hw+step1);
          float x3 = xy3*cos(gw);
          float y3 = xy3*sin(gw);
          float x4 = xy3*cos(gw+step2);
          float y4 = xy3*sin(gw+step2);

          X_normal_and_Vertex(groesse,x1,y1,z1);
          X_normal_and_Vertex(groesse,x2,y2,z1);
          X_normal_and_Vertex(groesse,x3,y3,z3);

          X_normal_and_Vertex(groesse,x2,y2,z1);
          X_normal_and_Vertex(groesse,x3,y3,z3);
          X_normal_and_Vertex(groesse,x4,y4,z3);


        
        }
      }
        glEnd();											// Done drawing the pyramid
      glDisable(GL_COLOR_MATERIAL);
      }


			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			System::Void OpenGL_Base::Render_House(System::Void)
			{
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
				glLoadIdentity();

				dd x, y, z, angle, l, b, h, x1, y1;
				x = 0.0;
				y = 0.0;
				z = 0.0;
				angle = 0.0;
				l = 3.0;
				b = 0.36;
				h = 2.0;

				Basic_Box(x, y, z, angle, l, b, h, x1, y1);


				x = x1;
				y = y1;
				l = 2.0;
				angle = 90.0;
				Basic_Box(x, y, z, angle, 2.0, b, h, x1, y1);
				x = x1;
				y = y1;
				Basic_Box_Fenster(x, y, z, angle, 1.0, b, h, 1.1, 1.8, x1, y1);
				x = x1;
				y = y1;
				Basic_Box(x, y, z, angle, 1.0, b, h, x1, y1);

				x = x1;
				y = y1;
				l = 1.0;
				angle = 45.0;
				Basic_Box(x, y, z, angle, l, b, h, x1, y1);

				x = x1;
				y = y1;
				l = 2.0;
				angle = 135.0;
				Basic_Box(x, y, z, angle, l, b, h, x1, y1);

				x = 5;
				y = 5;
				l = 10.0;
				Basic_Box(x, y, 0, 0, l, b, h, x1, y1); x = x1; y = y1;
				Basic_Box(x, y, 0, 90, l, b, h, x1, y1); x = x1; y = y1;
				Basic_Box(x, y, 0, 180, l, b, h, x1, y1); x = x1; y = y1;
				Basic_Box(x, y, 0, 270, l, b, h, x1, y1); x = x1; y = y1;


			}
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			System::Void OpenGL_Base::Render_House1(System::Void)
			{
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
				glLoadIdentity();

				dd x, y, z, angle, l, b, h, x1, y1;
				x = 0.0;
				y = 0.0;
				z = 0.0;
				angle = 0.0;
				l = 3.0;
				b = 0.36;
				h = 2.0;

				Basic_Box(x, y, z, angle, l, b, h, x1, y1);


				x = x1;
				y = y1;
				l = 2.0;
				angle = 90.0;
				Basic_Box(x, y, z, angle, 2.0, b, h, x1, y1);
				x = x1;
				y = y1;
				Basic_Box_Fenster(x, y, z, angle, 1.0, b, h, 1.1, 1.8, x1, y1);
				x = x1;
				y = y1;
				Basic_Box(x, y, z, angle, 1.0, b, h, x1, y1);

				x = x1;
				y = y1;
				l = 1.0;
				angle = 45.0;
				Basic_Box(x, y, z, angle, l, b, h, x1, y1);

				x = x1;
				y = y1;
				l = 2.0;
				angle = 135.0;
				Basic_Box(x, y, z, angle, l, b, h, x1, y1);

				x = 5;
				y = 5;
				l = 10.0;
				Basic_Box(x, y, 0, 0, l, b, h, x1, y1); x = x1; y = y1;
				Basic_Box(x, y, 0, 90, l, b, h, x1, y1); x = x1; y = y1;
				Basic_Box(x, y, 0, 180, l, b, h, x1, y1); x = x1; y = y1;
				Basic_Box(x, y, 0, 270, l, b, h, x1, y1); x = x1; y = y1;


			}
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			void OpenGL_Base::color1()
			{
				glColor3d(1.0, 1.0, 1.0);
			}
			void OpenGL_Base::color2()
			{
				glColor3d(0.9, 0.2, 0.2);
			}
			////////////////////////////////////////////////////////////////////////////////////////////
      System::Void OpenGL_Base::Basic_Box(dd x, dd y, dd z, dd angle, dd l, dd b, dd h, dd& x1, dd& y1)
      {
        glLoadIdentity();									// Reset the current modelview matrix
        glTranslated(x, y, z);
        glRotated(angle, 0, 0, 1);
        glBegin(GL_QUADS);									// Draw a quad
				color1();
				glVertex3d(0, 0, 0);                // Boden
        glVertex3d(0, b, 0);
        glVertex3d(l, b, 0);
        glVertex3d(l, 0, 0);

				color1();
				glVertex3d(0, 0, 0);         //Vorne
				glVertex3d(l, 0, 0);
				color2();
        glVertex3d(l, 0, h);
				color1();
        glVertex3d(0, 0, h);

				color1();
				glVertex3d(0, b, 0);         //Hinten
        glVertex3d(0, b, h);
				color2();
        glVertex3d(l, b, h);
        glVertex3d(l, b, 0);
				color1();

				color1();
				glVertex3d(0, 0, h);                //Oben
				color2();
        glVertex3d(l, 0, h);
        glVertex3d(l, b, h);
				color1();
        glVertex3d(0, b, h);

				color1();
        glVertex3d(0, 0, 0);                // Rechts
        glVertex3d(0, 0, h);
        glVertex3d(0, b, h);
        glVertex3d(0, b, 0);

				color1();
        glVertex3d(l, 0, 0);                // Links
        glVertex3d(l, b, 0);
				color2();
        glVertex3d(l, b, h);
        glVertex3d(l, 0, h);


        glEnd();//QUADS

        dd ar = 0.01745329251994329576923690768489 * angle;

        x1 = x + l*cos(ar);
        y1 = y + l*sin(ar);

      }
      //
      System::Void OpenGL_Base::Basic_Box_Fenster(dd x, dd y, dd z, dd angle, dd l, dd b, dd h, dd h1, dd h2, dd& x1, dd& y1)
      {
        Basic_Box(x, y, z, angle, l, b, h1, x1, y1);
        dd xx, yy;
        Basic_Box(x, y, h2, angle, l, b, h-h2, xx, yy);
      }
      //--------------------------------------------------------------
      System::Void OpenGL_Base::Basic_Disk(dd x, dd y, dd r)
      {

      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////
			System::Void OpenGL_Base::Render_Cube(System::Void)
			{
				{
					glDisable(GL_LIGHTING);

					glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
					glLoadIdentity();									// Reset the current modelview matrix

					glTranslatef(-1.5f, 0.0f, -0.0f);						// Move left 1.5 units and into the screen 6.0
					glRotatef(rtri_angle, 0.0f, 1.0f, 0.0f);						// Rotate the triangle on the y axis
					glBegin(GL_TRIANGLES);								// Start drawing a triangle
					glColor3f(1.0f, 0.0f, 0.0f);						// Red
					glVertex3f(0.0f, 1.0f, 0.0f);					// Top Of triangle (front)
					glColor3f(0.0f, 1.0f, 0.0f);						// Green
					glVertex3f(-1.0f, -1.0f, 1.0f);					// Left of triangle (front)
					glColor3f(0.0f, 0.0f, 1.0f);						// Blue
					glVertex3f(1.0f, -1.0f, 1.0f);					// Right of triangle (front)
					glColor3f(1.0f, 0.0f, 0.0f);						// Red
					glVertex3f(0.0f, 1.0f, 0.0f);					// Top Of triangle (right)
					glColor3f(0.0f, 0.0f, 1.0f);						// Blue
					glVertex3f(1.0f, -1.0f, 1.0f);					// Left of triangle (right)
					glColor3f(0.0f, 1.0f, 0.0f);						// Green
					glVertex3f(1.0f, -1.0f, -1.0f);					// Right of triangle (right)
					glColor3f(1.0f, 0.0f, 0.0f);						// Red
					glVertex3f(0.0f, 1.0f, 0.0f);					// Top Of triangle (back)
					glColor3f(0.0f, 1.0f, 0.0f);						// Green
					glVertex3f(1.0f, -1.0f, -1.0f);					// Left of triangle (back)
					glColor3f(0.0f, 0.0f, 1.0f);						// Blue
					glVertex3f(-1.0f, -1.0f, -1.0f);					// Right of triangle (back)
					glColor3f(1.0f, 0.0f, 0.0f);						// Red
					glVertex3f(0.0f, 1.0f, 0.0f);					// Top Of triangle (left)
					glColor3f(0.0f, 0.0f, 1.0f);						// Blue
					glVertex3f(-1.0f, -1.0f, -1.0f);					// Left of triangle (left)
					glColor3f(0.0f, 1.0f, 0.0f);						// Green
					glVertex3f(-1.0f, -1.0f, 1.0f);					// Right of triangle (left)
					glEnd();											// Done drawing the pyramid

					glLoadIdentity();									// Reset the current modelview matrix
					glTranslatef(1.5f, 0.0f, -0.0f);						// Move right 1.5 units and into the screen 7.0
					glRotatef(rquad_angle, 1.0f, 1.0f, 1.0f);					// Rotate the quad on the x axis 
					glBegin(GL_QUADS);									// Draw a quad
					glColor3f(0.0f, 1.0f, 0.0f);						// Set The color to green
					glVertex3f(1.0f, 1.0f, -1.0f);					// Top Right of the quad (top)
					glVertex3f(-1.0f, 1.0f, -1.0f);					// Top Left of the quad (top)
					glVertex3f(-1.0f, 1.0f, 1.0f);					// Bottom left of the quad (top)
					glVertex3f(1.0f, 1.0f, 1.0f);					// Bottom right of the quad (top)
					glColor3f(1.0f, 0.5f, 0.0f);						// Set The color to orange
					glVertex3f(1.0f, -1.0f, 1.0f);					// Top Right of the quad (bottom)
					glVertex3f(-1.0f, -1.0f, 1.0f);					// Top Left of the quad (bottom)
					glVertex3f(-1.0f, -1.0f, -1.0f);					// Bottom left of the quad (bottom)
					glVertex3f(1.0f, -1.0f, -1.0f);					// Bottom right of the quad (bottom)
					glColor3f(1.0f, 0.0f, 0.0f);						// Set The color to red
					glVertex3f(1.0f, 1.0f, 1.0f);					// Top Right of the quad (front)
					glVertex3f(-1.0f, 1.0f, 1.0f);					// Top Left of the quad (front)
					glVertex3f(-1.0f, -1.0f, 1.0f);					// Bottom left of the quad (front)
					glVertex3f(1.0f, -1.0f, 1.0f);					// Bottom right of the quad (front)
					glColor3f(1.0f, 1.0f, 0.0f);						// Set The color to yellow
					glVertex3f(1.0f, -1.0f, -1.0f);					// Top Right of the quad (back)
					glVertex3f(-1.0f, -1.0f, -1.0f);					// Top Left of the quad (back)
					glVertex3f(-1.0f, 1.0f, -1.0f);					// Bottom left of the quad (back)
					glVertex3f(1.0f, 1.0f, -1.0f);					// Bottom right of the quad (back)
					glColor3f(0.0f, 0.0f, 1.0f);						// Set The color to blue
					glVertex3f(-1.0f, 1.0f, 1.0f);					// Top Right of the quad (left)
					glVertex3f(-1.0f, 1.0f, -1.0f);					// Top Left of the quad (left)
					glVertex3f(-1.0f, -1.0f, -1.0f);					// Bottom left of the quad (left)
					glVertex3f(-1.0f, -1.0f, 1.0f);					// Bottom right of the quad (left)
					glColor3f(1.0f, 0.0f, 1.0f);						// Set The color to violet
					glVertex3f(1.0f, 1.0f, -1.0f);					// Top Right of the quad (right)
					glVertex3f(1.0f, 1.0f, 1.0f);					// Top Left of the quad (right)
					glVertex3f(1.0f, -1.0f, 1.0f);					// Bottom left of the quad (right)
					glVertex3f(1.0f, -1.0f, -1.0f);					// Bottom right of the quad (right)
					glEnd();											// Done drawing the quad
				}

			}
			////////////////////////////////////////////////////////////////////////////////////////////
			void X_Vertex(v3 c)
			{
				glVertex3f(c.x, c.y, c.z);
			}
			////////////////////////////////////////////////////////////////////////////////////////////
			void X_Vertex(cv3 c)
			{
				glVertex3f(c.x, c.y, c.z);
			}

			void X_Normal_Vertex(cv3 n, cv3 c)
			{
				glNormal3d(n.x,n.y,n.z);
				glVertex3d(c.x, c.y, c.z);
			}
			////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////
			System::Void render_arc(cv3& p0, cv3& v, cv3& q, double grad_v, double grad_q, double rmax, double rmin)
			{
				double rad_q = 3.1415*grad_q / 180.0;
				double rad_v = 3.1415*grad_v / 180.0;
				cv3 r = q%v; //richtung zum mittelpunkt
				r = cv3::drehung_vektor(r, v, rad_v);
				q = v % r;
				r.Normieren(1.0);
				q.Normieren(1.0);
				cv3 m0 = p0 + rmax*r;

				int N1 = 25;
				int N2 = 20;
				double delta_q = rad_q / N1;

				if (false)
				{
					glBegin(GL_TRIANGLES);								// Start drawing a triangle

					for (int i = 0; i < N1; i++)
					{
						cv3 p1 = cv3::drehung_punkt(p0, q, delta_q*i, m0);
						cv3 p2 = cv3::drehung_punkt(p0, q, delta_q*(i + 1), m0);
						glColor3f(1.0f, 0.0f, 0.0f);						// Red
						X_Vertex(p1);
						glColor3f(0.0f, 1.0f, 0.0f);						// Green
						X_Vertex(p2);
						glColor3f(0.0f, 0.0f, 1.0f);						// Blue
						X_Vertex(m0);
					}
					glEnd();//TRIANGLES
				}


				////////////////////////////////////////////////////////////////////
				glBegin(GL_TRIANGLES);								// Start drawing the tube
				for (int i = 0; i < N1; i++)
				{
					cv3 p1 = cv3::drehung_punkt(p0, q, delta_q*i, m0);
					cv3 p2 = cv3::drehung_punkt(p0, q, delta_q*(i + 1), m0);

					cv3 q1 = m0 - p1; q1.Normieren(1.0);
					cv3 q2 = m0 - p2; q1.Normieren(1.0);

					double delta_rmin = 6.28 / N2;
					for (int j = 0; j < N2; j++)
					{
						cv3 na = cos(delta_rmin*j)*q + sin(delta_rmin*j)*q1;
						cv3 nb = (cos(delta_rmin*(j + 1)))*q + (sin(delta_rmin*(j + 1)))*q1;
						cv3 nc = (cos(delta_rmin*(j + 1)))*q + (sin(delta_rmin*(j + 1)))*q2;
						cv3 nd = (cos(delta_rmin*j))*q + (sin(delta_rmin*j))*q2;

						cv3 a = p1 + rmin*na;
						cv3 b = p1 + rmin*nb;
						cv3 c = p2 + rmin*nc;
						cv3 d = p2 + rmin*nd;

						if (j % 2 == 0)
							glColor3f(1.0f, 0.0f, 0.0f);
						else
							glColor3f(0.0f, 1.0f, 0.0f);						// Green
						X_Normal_Vertex(na, a);
						X_Normal_Vertex(nb, b);
						X_Normal_Vertex(nc, c);
						//glColor3f(0.0f, 1.0f, 0.0f);						// Green
						X_Normal_Vertex(na, a);
						X_Normal_Vertex(nc, c);
						X_Normal_Vertex(nd, d);

					}


				}
				glEnd();//TRIANGLES

				p0 = cv3::drehung_punkt(p0, q, rad_q, m0);
				v = cv3::drehung_vektor(v, q, rad_q);
			}
			////////////////////////////////////////////////////////////////////////////////////////////
			System::Void render_cylinder(cv3& p0, cv3& v, cv3& q, double grad_v, double rmax, double rmin)
			{
				double rad_v = 3.1415*grad_v / 180.0;
				cv3 r = q%v; //richtung zum mittelpunkt
				r = cv3::drehung_vektor(r, v, rad_v);
				q = v % r;
				r.Normieren(1.0);
				q.Normieren(1.0);
				v.Normieren(1.0);

				int N1 = 20;
				int N2 = 15;
				double delta_qq = rmax / N1;
				////////////////////////////////////////////////////////////////////
				glBegin(GL_TRIANGLES);								// Start drawing the tube
				for (int i = 0; i < N1; i++)
				{
					cv3 p1 = p0 + (delta_qq*i)*v;
					cv3 p2 = p0 + (delta_qq*(i+1))*v;

					cv3 q1 = r;
					cv3 q2 = r;

					double delta_rmin = 6.28 / N2;
					for (int j = 0; j < N2; j++)
					{
						cv3 na = cos(delta_rmin*j)*q + sin(delta_rmin*j)*q1;
						cv3 nb = (cos(delta_rmin*(j + 1)))*q + (sin(delta_rmin*(j + 1)))*q1;
						cv3 nc = (cos(delta_rmin*(j + 1)))*q + (sin(delta_rmin*(j + 1)))*q2;
						cv3 nd = (cos(delta_rmin*j))*q + (sin(delta_rmin*j))*q2;

						cv3 a = p1 + rmin*na;
						cv3 b = p1 + rmin*nb;
						cv3 c = p2 + rmin*nc;
						cv3 d = p2 + rmin*nd;

						X_Normal_Vertex(na, a);
						X_Normal_Vertex(nb, b);
						X_Normal_Vertex(nc, c);
						//glColor3f(0.0f, 1.0f, 0.0f);						// Green
						X_Normal_Vertex(na, a);
						X_Normal_Vertex(nc, c);
						X_Normal_Vertex(nd, d);

					}


				}
				glEnd();//TRIANGLES

				p0 = p0 + (rmax*v);
			}


			
			////////////////////////////////////////////////////////////////////////////////////////////
			System::Void render_arc2(cv3& p0, cv3& v, cv3& q, double grad_v, double grad_q, double rmax, double rmin)
			{
				if (fabs(grad_q) < 0.01)
				{
					render_cylinder(p0, v, q, grad_v, rmax, rmin);
					return;
				}
				double rad_q = 3.1415*grad_q / 180.0;
				double rad_v = 3.1415*grad_v / 180.0;
				cv3 r = q%v; //richtung zum mittelpunkt
				r = cv3::drehung_vektor(r, v, rad_v);
				q = v % r;
				r.Normieren(1.0);
				q.Normieren(1.0);
				cv3 m0 = p0 + rmax*r;

				int N1 = 10;
				int N2 = 10;
				double delta_q = rad_q / N1;

				if (false)
				{
					glBegin(GL_TRIANGLES);								// Start drawing a triangle

					for (int i = 0; i < N1; i++)
					{
						cv3 p1 = cv3::drehung_punkt(p0, q, delta_q*i, m0);
						cv3 p2 = cv3::drehung_punkt(p0, q, delta_q*(i + 1), m0);
						glColor3f(1.0f, 0.0f, 0.0f);						// Red
						X_Vertex(p1);
						glColor3f(0.0f, 1.0f, 0.0f);						// Green
						X_Vertex(p2);
						glColor3f(0.0f, 0.0f, 1.0f);						// Blue
						X_Vertex(m0);
					}
					glEnd();//TRIANGLES
				}


				////////////////////////////////////////////////////////////////////
				glBegin(GL_TRIANGLES);								// Start drawing the tube
				for (int i = 0; i < N1; i++)
				{
					cv3 p1 = cv3::drehung_punkt(p0, q, delta_q*i, m0);
					cv3 p2 = cv3::drehung_punkt(p0, q, delta_q*(i + 1), m0);

					cv3 q1 = m0 - p1; q1.Normieren(1.0);
					cv3 q2 = m0 - p2; q2.Normieren(1.0);

					double delta_rmin = 6.28 / N2;
					for (int j = 0; j < N2; j++)
					{
						cv3 na = cos(delta_rmin*j)*q + sin(delta_rmin*j)*q1;
						cv3 nb = (cos(delta_rmin*(j + 1)))*q + (sin(delta_rmin*(j + 1)))*q1;
						cv3 nc = (cos(delta_rmin*(j + 1)))*q + (sin(delta_rmin*(j + 1)))*q2;
						cv3 nd = (cos(delta_rmin*j))*q + (sin(delta_rmin*j))*q2;

						cv3 a = p1 + rmin*na;
						cv3 b = p1 + rmin*nb;
						cv3 c = p2 + rmin*nc;
						cv3 d = p2 + rmin*nd;

						if (true)
						{
							if (j % 2 == 0)
								glColor3f(1.0f, 0.0f, 0.0f);
							else
								glColor3f(0.0f, 1.0f, 0.0f);						// Green
						}
						X_Normal_Vertex(na, a);
						X_Normal_Vertex(nb, b);
						X_Normal_Vertex(nc, c);
						//glColor3f(0.0f, 1.0f, 0.0f);						// Green
						X_Normal_Vertex(na, a);
						X_Normal_Vertex(nc, c);
						X_Normal_Vertex(nd, d);

					}


				}
				glEnd();//TRIANGLES

				p0 = cv3::drehung_punkt(p0, q, rad_q, m0);
				v = cv3::drehung_vektor(v, q, rad_q);
			}
			////////////////////////////////////////////////////////////////////////////////////////////
			GLfloat RubberAmbient1[] = { 0.05, 0.0, 0.0 }; //rubber
			GLfloat RubberAmbient[] = { 0.25, 0.1, 0.1 }; //rubber
			GLfloat RubberDiffuse[] = { 0.5, 0.4, 0.4 }; //rubber
			GLfloat RubberSpecular[] = { 0.7, 0.04, 0.04 }; //rubber


			GLfloat YellowRubberAmbient[] = { 0.05, 0.05, 0.0 }; //rubber
			GLfloat YellowRubberAmbient1[] = { 0.25, 0.25, 0.1 }; //rubber
			GLfloat YellowRubberDiffuse[] = { 0.5, 0.5, 0.4 }; //rubber
			GLfloat YellowRubberSpecular[] = { 0.7, 0.7, 0.04 }; //rubber

			void Make_Material(int iNr)
			{
				GLfloat shine = 128.0*.078125;
				switch (iNr)
				{
				case 0:
				default:
					glMaterialfv(GL_FRONT, GL_AMBIENT, RubberAmbient);
					glMaterialfv(GL_FRONT, GL_DIFFUSE, RubberDiffuse);
					glMaterialfv(GL_FRONT, GL_SPECULAR, RubberSpecular);
					glMaterialfv(GL_FRONT, GL_SHININESS, &shine);
					break;
				case 1:
					glMaterialfv(GL_FRONT, GL_AMBIENT, YellowRubberAmbient);
					glMaterialfv(GL_FRONT, GL_DIFFUSE, YellowRubberDiffuse);
					glMaterialfv(GL_FRONT, GL_SPECULAR, YellowRubberSpecular);
					glMaterialfv(GL_FRONT, GL_SHININESS, &shine);
					break;

				}
			}


			double timer = 0.0;
			////////////////////////////////////////////////////////////////////////////////////////////
			System::Void render_arc3(cv3& p0, cv3& v, cv3& q, double grad_v, double grad_q, double rmax, double rminA, double rminB)
			{
				//timer += 0.01;
				double wert = sin(timer);
				if (fabs(grad_q) < 0.01)
				{
					double rmin = 0.5*(rminA + rminB);
					render_cylinder(p0, v, q, grad_v, rmax, rmin);
					return;
				}
				//grad_q = wert*grad_q;
				double rad_q = 3.1415*grad_q / 180.0;
				double rad_v = 3.1415*grad_v / 180.0;
				cv3 r = q%v; //richtung zum mittelpunkt
				r = cv3::drehung_vektor(r, v, rad_v);
				q = v % r;
				r.Normieren(1.0);
				q.Normieren(1.0);
				cv3 m0 = p0 + rmax*r;

				int N1 = 15;
				int N2 = 20;
				double delta_q = rad_q / N1;

				int iNr = 0;
				////////////////////////////////////////////////////////////////////
				glBegin(GL_TRIANGLES);								// Start drawing the tube
				for (int i = 0; i < N1; i++)
				{
					iNr = i % 2;
					cv3 p1 = cv3::drehung_punkt(p0, q, delta_q*i, m0);
					cv3 p2 = cv3::drehung_punkt(p0, q, delta_q*(i + 1), m0);

					cv3 q1 = m0 - p1; q1.Normieren(1.0);
					cv3 q2 = m0 - p2; q2.Normieren(1.0);

					double delta_rmin = 6.28 / N2;
					for (int j = 0; j < N2; j++)
					{
						cv3 na = cos(delta_rmin*j)*q + sin(delta_rmin*j)*q1;
						cv3 nb = (cos(delta_rmin*(j + 1)))*q + (sin(delta_rmin*(j + 1)))*q1;
						cv3 nc = (cos(delta_rmin*(j + 1)))*q + (sin(delta_rmin*(j + 1)))*q2;
						cv3 nd = (cos(delta_rmin*j))*q + (sin(delta_rmin*j))*q2;


						double rmin1 = (rminB - rminA)*i / double(N1)+rminA;
						double rmin2 = (rminB - rminA)*(i+1) / double(N1)+rminA;
						cv3 a = p1 + rmin1*na;
						cv3 b = p1 + rmin1*nb;
						cv3 c = p2 + rmin2*nc;
						cv3 d = p2 + rmin2*nd;

						if (iNr == 1)
							iNr = 2;
						else
							iNr = 1;

						//Make_Material(iNr);

						X_Normal_Vertex(na, a);
						X_Normal_Vertex(nb, b);
						X_Normal_Vertex(nc, c);
						//glColor3f(0.0f, 1.0f, 0.0f);						// Green
						X_Normal_Vertex(na, a);
						X_Normal_Vertex(nc, c);
						X_Normal_Vertex(nd, d);

					}


				}
				glEnd();//TRIANGLES

				p0 = cv3::drehung_punkt(p0, q, rad_q, m0);
				v = cv3::drehung_vektor(v, q, rad_q);
			}
			////////////////////////////////////////////////////////////////////////////////////////////

			GLfloat redDiffuseMaterial[] = { 1.0, 0.0, 0.0 }; //set the 					material to red
			GLfloat whiteSpecularMaterial[] = { 1.0, 1.0, 1.0 }; //set					the material to white
			GLfloat greenEmissiveMaterial[] = { 0.0, 1.0, 0.0 }; //set					the material to green
			GLfloat whiteSpecularLight[] = { 1.0, 1.0, 1.0 }; //set the 					light specular to white
			GLfloat blackAmbientLight[] = { 0.0, 0.0, 0.0 }; //set the 				light ambient to black
			GLfloat whiteDiffuseLight[] = { 1.0, 1.0, 1.0 }; //set the 					diffuse light to white
			GLfloat blankMaterial[] = { 0.0, 0.0, 0.0 }; //set the diffuse					light to white
			GLfloat mShininess[] = { 128 }; //set the shininess of the 					material

			GLfloat whiteAmbientLight[] = { 1.0, 1.0, 1.0 }; //set the 				light ambient to white


			void init(void) {
				glEnable(GL_DEPTH_TEST);
				glEnable(GL_LIGHTING);
				glEnable(GL_LIGHT0);
				glEnable(GL_LIGHT1);
				//glEnable(GL_COLOR_MATERIAL);
				glEnable(GL_BLEND);
				glEnable(GL_POLYGON_SMOOTH);
			}

			void light(void) 
			{
				glPushMatrix();
				glLoadIdentity();
				GLfloat lightPosition0[] = { 2.0f, 3.0f, 2.0f }; //light
				glLightfv(GL_LIGHT0, GL_POSITION, lightPosition0);
				glLightfv(GL_LIGHT0, GL_SPECULAR, whiteSpecularLight);
				glLightfv(GL_LIGHT0, GL_AMBIENT, whiteAmbientLight);
				glLightfv(GL_LIGHT0, GL_DIFFUSE, whiteDiffuseLight);

				GLfloat lightPosition1[] = { -2.0f, -3.0f, -2.0f }; //light
				glLightfv(GL_LIGHT1, GL_POSITION, lightPosition0);
				glLightfv(GL_LIGHT1, GL_SPECULAR, whiteSpecularLight);
				glLightfv(GL_LIGHT1, GL_AMBIENT, whiteAmbientLight);
				glLightfv(GL_LIGHT1, GL_DIFFUSE, whiteDiffuseLight);
				glPopMatrix();
			}

			cv3 Random_on_sphere(Random^ r)
			{
				double lam = r->NextDouble(); 
				lam = lam * 6.283185307179586476925286766559;//; - 3.1415926535897932384626433832795;
				double phi = r->NextDouble();
				phi = asin(2.0*phi - 1.0);
				cv3 w;
				w.z = sin(phi);
				w.x = cos(phi)*cos(lam);
				w.y = cos(phi)*sin(lam);
				return w;
			}



			static GLuint  baseList1 = 0, baseList2 = 0;

			System::Void OpenGL_Base::Render_Tree(System::Void)
			{

				{
					init();
					glClearColor(0.2, 0.4, 0.6, 1.0);
					glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
					glLoadIdentity();
					light();



					Make_Material(0);

					int seeed = 41115;

					//glEnable(GL_COLOR_MATERIAL);
					if ((iRender1or2==1)&&(!glIsList(baseList1)))
					{
						baseList1 = glGenLists(1);
						glNewList(baseList1, GL_COMPILE);
						{
							double rmin_fak = 0.8;
							Random^ fixRand = gcnew Random(seeed);
							for (int i = 0; i < 35; i++)
							{
								double rmin = 0.5;
								cv3 p0;
								cv3 v = Random_on_sphere(fixRand);
								cv3 q1 = Random_on_sphere(fixRand);
								cv3 q = v%q1; q.Normieren(1.0);
								for (int j = 1; j <= 12; j++)
								{
									double rmax = fixRand->NextDouble()*0.5 + 1.5;
									rmax = 1.0;
									render_arc3(p0, v, q, fixRand->NextDouble() * 180.0 - 90.0, fixRand->NextDouble() *90.0+30.0, rmax, rmin, rmin_fak*rmin);
									rmin = rmin_fak*rmin;
								}
							}
						}
						glEndList();
					}

					if ((iRender1or2 == 2) && (!glIsList(baseList2)))
					{
						baseList2 = glGenLists(1);
						glNewList(baseList2, GL_COMPILE);
						{
							double rmin_fak = 0.8;
							Random^ fixRand = gcnew Random(seeed);
							for (int i = 0; i < 35; i++)
							{
								double rmin = 0.25;
								cv3 p0;
								cv3 v = Random_on_sphere(fixRand);
								cv3 q1 = Random_on_sphere(fixRand);
								cv3 q = v%q1; q.Normieren(1.0);
								for (int j = 1; j <= 12; j++)
								{
									double rmax = fixRand->NextDouble()*0.5 + 1.5;
									rmax = 1.0;
									render_arc3(p0, v, q, fixRand->NextDouble() * 180.0 - 90.0, fixRand->NextDouble() *90.0 + 30.0, rmax, rmin, rmin_fak*rmin);
									rmin = rmin_fak*rmin;
								}
							}
						}
						glEndList();
					}
					if (iRender1or2==1)
					  glCallList(baseList1);
					if (iRender1or2 == 2)
						glCallList(baseList2);				}
			}
			////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////
			void renderCylinder(float x1, float y1, float z1, float x2, float y2, float z2, float radius, int subdivisions, GLUquadricObj *quadric)
			{
				float vx = x2 - x1;
				float vy = y2 - y1;
				float vz = z2 - z1;

				//handle the degenerate case of z1 == z2 with an approximation
				if (vz == 0)
					vz = .0001;

				float v = sqrt(vx*vx + vy*vy + vz*vz);
				float ax = 57.2957795*acos(vz / v);
				if (vz < 0.0)
					ax = -ax;
				float rx = -vy*vz;
				float ry = vx*vz;
				glPushMatrix();

				//draw the cylinder body
				glTranslatef(x1, y1, z1);
				glRotatef(ax, rx, ry, 0.0);
				gluQuadricOrientation(quadric, GLU_OUTSIDE);
				gluCylinder(quadric, radius, radius, v, subdivisions, 1);

				//draw the first cap
				gluQuadricOrientation(quadric, GLU_INSIDE);
				gluDisk(quadric, 0.0, radius, subdivisions, 1);
				glTranslatef(0, 0, v);

				//draw the second cap
				gluQuadricOrientation(quadric, GLU_OUTSIDE);
				gluDisk(quadric, 0.0, radius, subdivisions, 1);
				glPopMatrix();
			}
			void renderCylinder_convenient(float x1, float y1, float z1, float x2, float y2, float z2, float radius, int subdivisions)
			{
				//the same quadric can be re-used for drawing many cylinders
				GLUquadricObj *quadric = gluNewQuadric();
				gluQuadricNormals(quadric, GLU_SMOOTH);
				renderCylinder(x1, y1, z1, x2, y2, z2, radius, subdivisions, quadric);
				gluDeleteQuadric(quadric);
			}
			////////////////////////////////////////////////////////////////////////////////////////////

			System::Void gitter(int iAnz, double dist)
			{
				GLUquadricObj *qobj;
				qobj = gluNewQuadric();

				for (int ix = -iAnz; ix <= iAnz; ix++)
					for (int iy = -iAnz; iy <= iAnz; iy++)
					{
						glPushMatrix();
						glTranslatef(ix*dist, iy*dist, 0);
						gluSphere(qobj, 0.125*dist, 10, 10);
						glPopMatrix();
					}


			}

			System::Void OpenGL_Base::Render_Cones(System::Void)
			{
				
				
					glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
					glLoadIdentity();									// Reset the current modelview matrix

					GLUquadricObj *qobj;
					qobj = gluNewQuadric();
					glColor3f(1.0f, 0.0f, 0.0f);						// Red
					gluQuadricDrawStyle(qobj, GLU_FILL); /* flat shaded */
					gluQuadricNormals(qobj, GLU_FLAT);
					
					glColor3f(0.0f, 0.0f, 1.0f);						// Blue-Z
					glLoadIdentity();									// Reset the current modelview matrix
					//glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
					gluCylinder(qobj, 0.25, 0.0, 2.0, 15, 5);

					glColor3f(0.0f, 1.0f, 0.0f);						// Green-Y
					glLoadIdentity();									// Reset the current modelview matrix
					glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
					gluCylinder(qobj, 0.25, 0.0, 2.0, 15, 5);

					glColor3f(1.0f, 0.0f, 0.0f);						// Red-X
					glLoadIdentity();									// Reset the current modelview matrix
					glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
					gluCylinder(qobj, 0.25, 0.0, 2.0, 15, 5);


					glColor3f(1.0f, 1.0f, 0.0f);						// Yellow
					glLoadIdentity();									// Reset the current modelview matrix

					renderCylinder(2, 0, 0, 0, 2, 0, 0.25, 15, qobj);

					glPushMatrix();
					glTranslatef(2, 0, 0);
					gluSphere(qobj, 0.25, 15, 15);
					glPopMatrix();

					glPushMatrix();
					glTranslatef(0, 2, 0);
					gluSphere(qobj, 0.25, 15, 15);
					glPopMatrix();


					double ll = 3.0;
					glLoadIdentity();
					glLineWidth(15.0f);
					glBegin(GL_LINES);
					glColor3f(0.5f, 0.0f, 0.00f);						// Red
					X_Vertex(ll, 0.0, 0.0, 0.0);
					X_Vertex(ll, 1.0, 0.0, 0.0);
					glColor3f(0.0f, 0.5f, 0.00f);						// g
					X_Vertex(ll, 0.0, 0.0, 0.0);
					X_Vertex(ll, 0.0, 1.0, 0.0);
					glColor3f(0.0f, 0.0f, 0.5f);						// b
					X_Vertex(ll, 0.0, 0.0, 0.0);
					X_Vertex(ll, 0.0, 0.0, 1.0);
					glEnd();

					//Kugeln
					glLoadIdentity();
					glColor3f(0.9f, 0.0f, 0.00f);						// Red
					gitter(10, 1.0);

					glPushMatrix();
					glRotatef(90.0, 1.0, 0.0, 0.0);
					glColor3f(0.0f, 0.9f, 0.00f);						// Red
					gitter(10, 1.0);
					glPopMatrix();

					glPushMatrix();
					glRotatef(90.0, 0.0, 1.0, 0.0);
					glColor3f(0.0f, 0.0f, 0.900f);						// Red
					gitter(10, 1.0);
					glPopMatrix();

			}

			double fkt(double gw, double hw)
			{
				double pihalbe = 0.5*3.1415926535;
				int iTyp = 0;
				switch (iTyp)
				{
				case 1:
					return hw*hw + gw*gw;
					break;
				case 2:
					return atan2(hw, gw);
					break;
				case 3:
					return atan2(gw, 3 * hw);
					break;
				case 4:
					return 3.0*(1.0 + 0.1*sin(5.0*(hw + gw) + 0.15*cos(6.0*(hw - gw)))+0.075*cos(5*hw*hw+4*hw*gw+3*gw*gw))/(1+hw*hw);
					break;
				case 5:
					return 4.0/(0.01+(hw-pihalbe)*(-hw-pihalbe));
					break;
				}
					return 1.5*(1.0 + 0.1*sin(5.0*(hw + gw)));
			}
			////////////////////////////////////////////////////////////////////////////////////////////
			System::Void OpenGL_Base::Render_WaveBall(System::Void)
			{
				double pihalbe = 0.5*3.1415926535;
				double dStep = pihalbe / 25.0;
				v3 a, b, c, d;

				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
				glLoadIdentity();									// Reset the current modelview matrix

				glBegin(GL_TRIANGLES);
				int ii = 0;
				for (double hw = -pihalbe+0.5; hw <= pihalbe-0.5; hw += dStep)
				{
					ii+=2;
					for (double gw = 0.0; gw <= 4.0*pihalbe; gw += dStep)
					{


						a.set_polar(gw, hw, fkt(gw,hw));
						b.set_polar(gw + dStep, hw, fkt(gw+dStep, hw));
						c.set_polar(gw + dStep, hw + dStep, fkt(gw + dStep, hw + dStep));
						d.set_polar(gw, hw + dStep, fkt(gw, hw + dStep));

						ii++;
						if (ii % 2 == 0)
							glColor3f(1.0f, 0.0f, 0.0f);						// Red
						else
							glColor3f(0.0f, 0.0f, 1.0f);						// Blue

						X_Vertex(a);
						X_Vertex(b);
						X_Vertex(c);

						X_Vertex(a);
						X_Vertex(c);
						X_Vertex(d);

					}
				}
				glEnd();
			}
			////////////////////////////////////////////////////////////////////////////////////////////
		System::Void OpenGL_Base::Render_Gelände(System::Void)
		{

			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
			glLoadIdentity();									// Reset the current modelview matrix
			glTranslatef(-1.5f, 0.0f, -0.0f);						// Move left 1.5 units and into the screen 6.0

			glBegin(GL_TRIANGLES);								// Start drawing a triangle

			float step = 0.25;
			float groesse = 10.0;
			int ii = 0;
			for (float x = -groesse; x <= groesse; x += step)
				for (float y = -groesse; y <= groesse; y += step)
				{
					ii++;
					if (ii % 2 == 0)
						glColor3f(1.0f, 0.0f, 0.0f);						// Red
					else
						glColor3f(0.0f, 0.0f, 1.0f);						// Blue

					//glColor3f(1.0f,0.0f,0.0f);						// Red
					glVertex3f(x, fkt(berg_angle, x, y), y);					// Top Of triangle (front)
					//glColor3f(1.0f,0.0f,0.0f);						// Red
					glVertex3f(x + step, fkt(berg_angle, x + step, y), y);						// Left of triangle (front)
					//glColor3f(1.0f,0.0f,0.0f);						// Red
					glVertex3f(x, fkt(berg_angle, x, y + step), y + step);						// Right of triangle (front)

					//glColor3f(0.0f,0.0f,1.0f);						// Blue
					glVertex3f(x + step, fkt(berg_angle, x + step, y), y);					// Top Of triangle (front)
					//glColor3f(0.0f,0.0f,1.0f);						// Blue
					glVertex3f(x + step, fkt(berg_angle, x + step, y + step), y + step);						// Left of triangle (front)
					//glColor3f(0.0f,0.0f,1.0f);						// Blue
					glVertex3f(x, fkt(berg_angle, x, y + step), y + step);						// Right of triangle (front)


				}
			glEnd();											// Done drawing the pyramid
		}
		////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////
		System::Void OpenGL_Base::Render_Material(System::Void)
		{

			// set up light
			float ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
			float diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
			float specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
			float position[] = { 200.0f, 300.0f, 100.0f, 1.0f };

			glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
			glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
			glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
			glLightfv(GL_LIGHT0, GL_POSITION, position);
			glEnable (GL_LIGHT0);
			glEnable (GL_LIGHTING);


			float mat_ambient1[] = { 0.129412f, 0.223529f, 0.027451f };
			float mat_diffuse1[] = { 0.780392f, 0.568627f, 0.113725f };
			float mat_specular1[] = { 0.992157f, 0.941176f, 0.807843f };
			float shine1 = 13.8f;

			//Emerald
			float mat_ambient2[] = { 0.0215f, 0.1745f, 0.9215f, 0.55f };
			float mat_ambient2a[] = { 0.0215f, 0.1745f, 0.0215f, 0.55f };
			float mat_diffuse2[] = { 0.07568f, 0.61424f, 0.07568f, 0.55f };
			float mat_specular2[] = { 0.633f, 0.727811f, 0.633f, 0.55f };
			float shine2 = 76.8f;

			glEnable(GL_COLOR_MATERIAL);
			glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);



			glShadeModel(GL_FLAT);




			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
			glLoadIdentity();									// Reset the current modelview matrix
			glTranslatef(-1.5f, 0.0f, -0.0f);						// Move left 1.5 units and into the screen 6.0

			glBegin(GL_TRIANGLES);								// Start drawing a triangle


      float nx,ny,nz;
      float h = 0.1;

			float step = 0.25;
			float groesse = 5.0;
			int ii = 0;
			for (float x = -groesse; x <= groesse; x += step)
				for (float y = -groesse; y <= groesse; y += step)
				{
					ii++;
					if (ii % 2 == 0)
					{
						
						glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse1);
						glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular1);
						glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient1);
						glMaterialf(GL_FRONT, GL_SHININESS, shine1);
						
				//glColor3f(0.0f,1.0f,0.0f);						// Green
						glEnable(GL_COLOR_MATERIAL);

					}
					else
					{
						glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
						glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse2);
						glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular2);
						glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient2);
						glMaterialf(GL_FRONT, GL_SHININESS, shine2);
						glEnable(GL_COLOR_MATERIAL);
					}

					//glColor3f(1.0f,0.0f,0.0f);						// Red
          normale_fkt(berg_angle, x, y, fkt,  h, nx,  ny,  nz);
          glNormal3f(nx,nz,ny);
					glVertex3f(x, fkt(berg_angle, x, y), y);					// Top Of triangle (front)
					//glColor3f(1.0f,0.0f,0.0f);						// Red
          normale_fkt(berg_angle, x+step, y, fkt,  h, nx,  ny,  nz);
          glNormal3f(nx,nz,ny);
					glVertex3f(x + step, fkt(berg_angle, x + step, y), y);						// Left of triangle (front)
					//glColor3f(1.0f,0.0f,0.0f);						// Red
          normale_fkt(berg_angle, x, y+step, fkt,  h, nx,  ny,  nz);
          glNormal3f(nx,nz,ny);
					glVertex3f(x, fkt(berg_angle, x, y + step), y + step);						// Right of triangle (front)

					//glColor3f(0.0f,0.0f,1.0f);						// Blue
          normale_fkt(berg_angle, x + step, y, fkt,  h, nx,  ny,  nz);
          glNormal3f(nx,nz,ny);
					glVertex3f(x + step, fkt(berg_angle, x + step, y), y);					// Top Of triangle (front)
					//glColor3f(0.0f,0.0f,1.0f);						// Blue
          normale_fkt(berg_angle, x + step, y + step, fkt,  h, nx,  ny,  nz);
          glNormal3f(nx,nz,ny);
					glVertex3f(x + step, fkt(berg_angle, x + step, y + step), y + step);						// Left of triangle (front)
					//glColor3f(0.0f,0.0f,1.0f);						// Blue
          normale_fkt(berg_angle, x, y + step, fkt,  h, nx,  ny,  nz);
          glNormal3f(nx,nz,ny);
					glVertex3f(x, fkt(berg_angle, x, y + step), y + step);						// Right of triangle (front)

					glDisable(GL_COLOR_MATERIAL);
				}
			glEnd();											// Done drawing the pyramid

			glDisable(GL_LIGHTING);
			glDisable(GL_LIGHT0);
			glDisable(GL_COLOR_MATERIAL);

		}
		////////////////////////////////////////////////////////////////////////////////////////////
		System::Void OpenGL_Base::Render_Material_A(System::Void)
		{
			//************** Object 1 **************
			glEnable(GL_COLOR_MATERIAL);
			glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
			glColor4f(149.0 / 255.0, 78.0 / 255.0, 22.0 / 255.0, 1.0);
			float mat_specular[] = { 0.992157, 0.941176, 0.807843, 1.0 };
			float shininess = 10;

			glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
			glMaterialf(GL_FRONT, GL_SHININESS, shininess);

			glPushMatrix();
			glTranslatef(0, 3.0, 0);
			//drawSphere(0.1, 0.1, 0.1);
			glRotatef(10, 1, 0, 0);

			glDisable(GL_COLOR_MATERIAL);


			//************** Object 2 *****************
			glEnable(GL_COLOR_MATERIAL);
			glColorMaterial(GL_FRONT, GL_DIFFUSE);
			glColor4f(48.0 / 255.0, 48.0 / 255.0, 48.0 / 255.0, 1.0);
			float mat_specular_2[] = { 0.992157, 0.941176, 0.807843, 1.0 };
			float shininess_2 = 10;

			glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular_2);
			glMaterialf(GL_FRONT, GL_SHININESS, shininess_2);

			glPushMatrix();
			glTranslatef(-0.6, 0.2, 1.6 / 2.0);
			//drawSphere(0.1, 0.1, 0.1);
			glPopMatrix();

			glDisable(GL_COLOR_MATERIAL);
		}
		////////////////////////////////////////////////////////////////////////////////////////////
		System::Void OpenGL_Base::Render(System::Void)
		{
      if(m_iSetPixelFormat1)
      {
      
      switch (m_iShowCubes_or_Gel)
      {
        case eCube:
					Render_Cube();
					//Render_Cones();
					break;
				case eTree:
					Render_Tree();
					//Render_Cones();
					break;
				case ePlaneten:
					//Render_House();
					Render_Planeten();
					//Render_Cones();
					break;
				case eCones:
					Render_Cones();
					//Render_Gelände();
					break;
				case eWaveSphere:
					Render_WaveBall();
					//Render_Gelände();
					break;
				case eMaterial:
					Render_Material();
					break;
				case eSphere_Material:
					Render_Sphere_Material();
					break;
				case eSphere:
					Render_Sphere();
					break;
				case eWurmloch:
					Render_Wurmloch();
					break;
				case eMinKugeln:
					Render_MinKugeln();
					break;
				case 72:
					Render_Archimedes();
					break;

					break;
			}
			rtri_angle+=0.2f;											// Increase the rotation variable for the triangle
			rquad_angle-=0.15f;										// Decrease the rotation variable for the quad
      berg_angle+= 0.01;
      }  
		}
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

		System::Void OpenGL_Base::SwapOpenGLBuffers1(System::Void)
		{
			if (m_iSetPixelFormat1)
			{
				SwapBuffers(m_HDC1);
			}
		}


		System::Void OpenGL_Base::SwapOpenGLBuffers2(System::Void)
		{
			if (m_iSetPixelFormat2)
			{
				SwapBuffers(m_HDC2);
			}
		}


		System::Void OpenGL_Base::ToggleShow(enumScenenAuswahl i)
   {
		 m_iShowCubes_or_Gel = i;
   }


	 System::Void OpenGL_Base::ReadObjFile(String^ path)
	 {
		 using namespace System::IO;

		 String^ line;
		 array<Char>^chars = { ' ', ',', '->', ':' };

		 m_Obj = new MyObj();

		 StreamReader^ sr = gcnew StreamReader(path);
		 try
		 {
			 while (sr->Peek() >= 0)
			 {
				 line = sr->ReadLine();
				 line->Trim();
				 //line->Replace("   ", " ");
				 //line->Replace("  ", " ");
				 //line->Replace("  ", " ");
				 line->Replace(".", ",");
				 if (line->Length > 0)
				 {
					 if (line->Substring(0, 1) == "v")
					 {
						 cli::array<String^>^ ar = line->Split(chars, StringSplitOptions::RemoveEmptyEntries);
						 if (ar->Length == 4 || ar->Length == 5)
						 {
							 double x = Convert::ToDouble(ar[1], System::Globalization::CultureInfo::InvariantCulture); //Convert.ToDouble(stringEnglish, CultureInfo.InvariantCulture);
							 double y = Convert::ToDouble(ar[2], System::Globalization::CultureInfo::InvariantCulture);
							 double z = Convert::ToDouble(ar[3], System::Globalization::CultureInfo::InvariantCulture);
							 m_Obj->m_v.add(x, y, z);
						 }
					 }
					 else if (line->Substring(0, 1) == "f")
					 {
						 cli::array<String^>^ ar = line->Split(chars, StringSplitOptions::RemoveEmptyEntries);
						 if (ar->Length >= 4)
						 {
							 int ifacecount = ar->Length - 1; ifacecount = -ifacecount;
							 m_Obj->m_f.push_back(ifacecount);
							 for (int i = 1; i < ar->Length; i++)
							 {
								 int iface = Convert::ToInt32(ar[i]);
								 m_Obj->m_f.push_back(iface);
							 }
						 }
					 }
				 }
			 }
		 }
		 finally
		 {
			 delete sr;
		 }

	 }
	 ////////////////////////////////////////////////////////////////////
	 ////////////////////////////////////////////////////////////////////
	 int MyObj::get_face(int ifaceNr) //nullbased
	 {
		 int ipos = 0;
		 int iface = 0;
		 while (iface < ifaceNr)
		 {
			 int iv = m_f[ipos];
			 if (iv < 0)
			 {
				 ipos += abs(iv);
				 iface++;
			 }
			 else
			 {
				 return -1;
			 }
		 }
		 return ipos;
	 }
	 //
	 cv3 MyObj::get_face_normal(int ifaceNr)
	 {
		 cv3 n;
		 int iface_pos = get_face(ifaceNr);
		 if (iface_pos >= 0)
		 {
			 int iv = m_f[iface_pos];
			 if (iv < 0)
			 {
				 cv3 a = m_v[m_f[iface_pos + 2]] - m_v[m_f[iface_pos + 1]];
				 cv3 b = m_v[m_f[iface_pos + 3]] - m_v[m_f[iface_pos + 1]];
				 n = a%b;
				 n.Normieren(1.0);
				 return n;
			 }
			 else
			 {
				 return n;
			 }
		 }

	 }


	 CGravitation::CGravitation()
	 {
		 p.add(1, 1, 1, 0.5);
		 p.add(1, -1, 0, 0.25);
		 p.add(0, 1, -1, 0.25);
		 v.add(0.01, 0, 0);
		 v.add(0., 0.01, 0);
		 v.add(0., 0, 0.01);
		 
		 f.add(0, 0, 0);
		 f.add(0, 0, 0);
		 f.add(0, 0, 0);

	 }

	 void CGravitation::Go()
	 {
		 for (int i = 0; i < f.size(); i++)
		 {
			 f[i].set_null();

			 cv3 fi;

		 }


		 for (int i = 0; i < v.size(); i++)
		 {
			 p[i] = p[i] + v[i];
		 }
	 }