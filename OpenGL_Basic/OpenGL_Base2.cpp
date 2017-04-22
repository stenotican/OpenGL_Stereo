#include "StdAfx.h"
#include "OpenGL_Base.h"



#include "nu_optimierung.h"
#include "MyData.h"

#include <math.h>


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
List<v3^>^ Wurmloch::GetLine1poly(int N, float r1, float w1, float r2, float w2)
{
	List<v3^>^arr = gcnew List < v3^ >;

	CDoubleVector vec_start(N);

	nu_FunctorN_Wurmloch_Polygon fn;
	fn.t1 = r1;
	fn.t2 = r2;
	fn.w1 = w1;
	fn.w2 = w2;
	fn.N = N;
	fn.NLines = 50;


	nu_Optimierung opt;
	int NDim = vec_start.size();
	double ftol = 0.0000001;
	double lambda = 0.01;
	CDoubleVector  vec_erg(N);

	if (opt.nu_amoeba(NDim, ftol, lambda, vec_start, fn, vec_erg))
	{
		for (int i = 0; i <= fn.NLines; i++)
		{
			float t = float(i) / float(fn.NLines);
			float r = fn.polygon(t, vec_erg, r1, r2);
			float w = w1 + t*(w2 - w1);
			arr->Add(gcnew v3(r, w));

		}
	}
	else
	{
		for (int i = 0; i <= fn.NLines; i++)
		{
			float t = float(i) / float(fn.NLines);
			float r = r1 + t*(r2 - r1);
			float w = w1 + t*(w2 - w1);
			arr->Add(gcnew v3(r, w));

		}
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
	if (iArt == 2)
	{
		w2 = w2 - 6.283185307179586476925286766559; //=2 *3.1415926535897932384626433832795
	}


	std::vector<double> va(N + 1);
	std::vector<double> vt(N + 1);

	std::vector<cv3> vp(N + 1);

	for (int i = 0; i <= N; i++)
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

			double dd = max(d1a, d2a);

			cv3 m = 0.5*(p1 + p2);
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
			while (ic > 0 && dd2 > dd)
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
List<v3^>^ Wurmloch::GetLine6_steps(int iArt, int N, float r1, float w1, float r2, float w2)
{
	// r=h/cos(alpha)
	if (iArt == 2)
	{
		w2 = w2 - 6.283185307179586476925286766559; //=2 *3.1415926535897932384626433832795
	}


	std::vector<double> va(N + 1);
	std::vector<double> vt(N + 1);

	std::vector<cv3> vp(N + 1);

	cv3 v1(r1, w1);
	cv3 v2(r2, w2);

	double tsenkrecht = (v1 - v2)*v1;
	double v21 = (v2 - v1)*(v2 - v1);
	tsenkrecht = tsenkrecht/v21;

	cv3 vsenkrecht = v1 + tsenkrecht*(v2 - v1);

	double wsenkrecht = atan2(vsenkrecht.y, vsenkrecht.x);
	double hsenkrecht = sqrt(vsenkrecht.y*vsenkrecht.y + vsenkrecht.x*vsenkrecht.x);



	for (int i = 0; i <= N; i++)
	{
		double t = float(i) / float(N);
		double r = r1 + t*(r2 - r1);
		double w = w1 + t*(w2 - w1);
		if (r1*r2>0.0)
		  r = hsenkrecht / cos(w - wsenkrecht);
		va[i] = w;
		vt[i] = r;
		cv3 p(r, w);
		vp[i] = p;
	}


	int icount = 1000;
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

			double dd = max(d1a, d2a);

			cv3 m = 0.5*(p1 + p2);
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
			while (ic > 0 && dd2 > dd)
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
List<v3^>^ Wurmloch::MinKugeln(int iArt, int N, int iAnzIter) //kugeln auf einheitskugel
{
	CDoubleVector vec_start(2*N-3);

	nu_FunctorN_MinKugel fn;
	fn.NN = N;
	fn.initVec(iArt,vec_start);
	fn.m_iOptiZiel = 0;
	nu_Optimierung opt;
	opt.NMAX = iAnzIter;
	int NDim = vec_start.size();
	double ftol = 1.0E-9;
	double lambda = 0.1;
	CDoubleVector  vec_erg(NDim);
	vec_erg = vec_start;

	int iok = 0;
	if (opt.nu_amoeba(NDim, ftol, lambda, vec_start, fn, vec_erg))
	{
		iok = 1;
	}

	List<v3^>^arr = gcnew List < v3^ >;

	fn.reset_vec_cv3(vec_erg);
	MinDist = fn.mindist();

	for (int i = 0; i < N; i++)
	{
		arr->Add(gcnew v3(fn.v_cv3[i].x, fn.v_cv3[i].y, fn.v_cv3[i].z));
	}
	return arr;

}
//
//
//
//
//
//
//
//
OpenGL_Base::OpenGL_Base(System::Windows::Forms::Panel^  panel1, System::Windows::Forms::Panel^  panel2)
{
	//
	//double a11=-2.0, a12=5.0, a21=3.0, a22=-4.0, b1=29.0, b2=-26.0, x1, x2;
	//solve2( a11,  a12,  a21,  a22,  b1,  b2,  x1,  x2);
  //
  m_iDoMouseMove = 0;
  m_iShowCubes_or_Gel=1;

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
///////////////////////////////////////////////////
bool OpenGL_Base::InitGL(GLvoid)										// All setup for opengl goes here
		{
			glShadeModel(GL_SMOOTH);							// Enable smooth shading
			glClearColor(0.5f, 0.5f, 0.7f, 0.5f);				// Black background
			glClearDepth(1.0f);									// Depth buffer setup
			glEnable(GL_DEPTH_TEST);							// Enables depth testing
			glDepthFunc(GL_LEQUAL);								// The type of depth testing to do
			glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really nice perspective calculations
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


	void X_Vertex(v3 c)
	{
		glVertex3f(c.x, c.y, c.z);
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
				int NPoly = 50;
				if (arr_glob1 == nullptr)
				{
					arr_glob1 = Wurmloch::GetLine6_steps(1,NPoly, 4.0f, 0.0f, 2.0f, 2.8f);
				}
				if (arr_glob2 == nullptr)
				{
					arr_glob2 = Wurmloch::GetLine6_steps(1,NPoly, 4.0f, 0.0f, -2.0f, 2.8f);
				}
				if (arr_glob3 == nullptr)
				{
					arr_glob3 = Wurmloch::GetLine6_steps(2,NPoly, 4.0f, 0.0f, -2.0f, 2.8f);
				}

				List<v3^>^ arr = arr_glob1;
				glColor3f(0.99f, 1.0f, 0.01f);						// Yellow
				glLineWidth(15.0f);
				glBegin(GL_LINES);
				int NN = arr->Count - 1;
				for (int i = 0; i < NN; i++)
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
				for (int i = 0; i < NN; i++)
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
				for (int i = 0; i < NN; i++)
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
			case 5:
			{
				int N = 75;
				if (arr_glob1 == nullptr)
				{
					arr_glob1 = Wurmloch::GetLine5_steps(1, N, 4.0f, 0.0f, 4.0f, 2.0f);
				}
				if (arr_glob2 == nullptr)
				{
					arr_glob2 = Wurmloch::GetLine5_steps(1, N, 4.0f, 0.0f, -4.0f, 2.0f);
				}
				if (arr_glob3 == nullptr)
				{
					arr_glob3 = Wurmloch::GetLine5_steps(2, N, 4.0f, 0.0f, -4.0f, 2.0f);
				}

				List<v3^>^ arr = arr_glob1;
				glColor3f(0.99f, 1.0f, 0.01f);						// Yellow
				glLineWidth(15.0f);
				glBegin(GL_LINES);
				for (int i = 0; i < N; i++)
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
				for (int i = 0; i < N; i++)
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
				for (int i = 0; i < N; i++)
				{
					v3^ a = arr[i];
					v3^ b = arr[i + 1];
					glVertex3f(a->x, a->y, a->z);
					glVertex3f(b->x, b->y, b->z);
				}
				glEnd();



			}
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
				arr_glob1 = Wurmloch::MinKugeln(iFormFaktor, N, m_data.m_iMaxIterations);
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
	////////////////////////////////////////////////////////////////////////////////////////////
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
			////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////
			System::Void OpenGL_Base::cmBox(float x, float y, float w, float dx, float dy, float dz)
			{
				float n0 = 0.0f;
				float n1 = 1.0f;
				glColor3f(0.0f, 1.0f, 0.0f);						// Set The color to green

				glBegin(GL_QUADS);									// Draw a box
				glVertex3f(n0, n0, n0);				
				glVertex3f(dx, n0, n0);
				glVertex3f(dx, dy, n0);
				glVertex3f(n0, dy, n0);
				glEnd();
			}
			////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////
			System::Void OpenGL_Base::Render_Haus(System::Void)
			{
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear screen and depth buffer
				glLoadIdentity();

				cmBox(2.0, 3.0, 45.0, 5.0, 0.5, 3.0);

			}


			System::Void OpenGL_Base::Render_Cube(System::Void)
			{
				{


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
        case 1:
					Render_Cube();
					//Render_Cones();
					break;
				case 2:
					Render_WaveBall();
					//Render_Gelände();
					break;
				case 3:
					//Render_Material();
					Render_Haus();
					break;
				case 4:
					Render_Sphere_Material();
					break;
				case 5:
					Render_Sphere();
					break;
				case 6:
					Render_Wurmloch();
					break;
				case 11:
					Render_MinKugeln();

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


   System::Void OpenGL_Base::ToggleShow(int i)
   {
		 m_iShowCubes_or_Gel = i;
   }

