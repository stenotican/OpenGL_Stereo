#pragma once

#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include <vector>

#include "TriangelTool.h"
#include "MyData.h"

#include "cv3.h"



using namespace System;
using namespace System::Windows::Forms;
using namespace System::Collections::Generic;



namespace OpenGLForm 
{

	enum enumScenenAuswahl{ eCube, eTree, eMinKugeln, ePlaneten, eCones, eSphere, eGelände, eWaveSphere, eSphere_Material, eMaterial, eWurmloch, eArchimedes };


	public ref class v3
	{
	public:
		float x, y, z;
	public:
		v3(){ x = 0.0; y = 0.0; z = 0.0; }
		v3(const v3% v){ x = v.x; y = v.y; z = v.z; }
		v3(double a, double b, double c){ x = a; y = b; z = c; }
		v3(float a, float b, float c){ x = a; y = b; z = c; }
		v3(float r, float w); //xyz-Wurm
		v3(int iArt,float r, float w); //xyz-Normale

		void set_wurm(float r, float w); //xyz-Wurm

		void set_polar(double gw, double hw, double r);

		void mal(double faktor){ x = faktor *x; y = faktor *y; z = faktor *z; }
		void set(double a, double b, double c){ x = a; y = b; z = c; }
		void set_p(v3^ vp);
		void set_r(v3% vr);

	public:
		float Dist(v3% v);
		float Norm();// { return sqrt(x*x + y*y + z*z); }
		void Normieren(float laenge);
		v3^ langs_this(const v3% a);
		v3^ quer_this(const v3% a);

		v3^ v3::operator+(const v3% a);
		v3^ v3::operator-(const v3% a);

		virtual System::String^ ToString() override;


	};

	//v3^ operator+(const v3% a, const v3% b);
	//v3^ operator-(const v3% a, const v3% b);
	static v3^ operator*(const float a, const v3% b);

	static float operator*(const v3% a, const v3% b);
	static v3^ operator%(const v3% a, const v3% b);
	//float operator*(const v3^ a, const v3^ b);

	//////////////////////////////////////////////////////////////
	public ref class i3
	{
	public:
		int ix, iy, iz;
		int% i(int i){ return ix; }
	public:
		i3(){ ix = iy = iz = 0; }
		i3(const i3% v){ ix = v.ix; iy = v.iy; iz = v.iz; }
		i3(int a, int b, int c){ ix = a; iy = b; iz = c; }
	};

	public ref class listv3 : System::Collections::Generic::List<v3^>
	{
	public:
		void add(String^ s3);
	};

	public ref class listi3 : System::Collections::Generic::List<i3^>
	{
	public:
		void add(String^ s3);
	};

	public ref class triangle
	{
	public:
		listv3 m_v; //vertex 
		listv3 m_n; //noraml
		listi3 m_e; //edge

	};

	///////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////



	public ref class Wurmloch
	{
	public:
		static double MinDist = 0.0;
		static int m_iOptiZiel = 0;

	public:
		static void get_rw(v3^v, float% r, float% w, float% dist, v3^% s3);
    static double PolynomWert(double x, std::vector<double>& ai);
    static void get_rw_polynom(double t01, double r1, double w1, double r2, double w2, std::vector<double>& ai, double& w, double& r);

	public: 
		static List<v3^>^ GetLine1(int N, float r1, float w1, float r2, float w2);
		//static List<v3^>^ GetLine1poly(int N, float r1, float w1, float r2, float w2);

		static List<v3^>^ GetLine2(int N, float r1, float w1, float r2, float w2);

		static List<v3^>^ GetLine3(int N, float r1, float w1, float r2, float w2);

		static List<v3^>^ MakeLine(std::vector<double>& arr_tw);
		static void       GetArr_TW(int N, float r1, float w1, float r2, float w2, std::vector<double>& arr_tw);

		static List<v3^>^ GetLine4Winkel(int N, float r1, float w1, float r2, float w2);
		static List<v3^>^ GetLine5_steps(int iArt, int N, float r1, float w1, float r2, float w2);
		//static List<v3^>^ GetLine6_steps(int iArt, int N, float r1, float w1, float r2, float w2);

		static List<v3^>^ MinKugeln(int iArt, int N, int iAnzIter, int iOptiZiel); //kugeln auf einheitskugel
	};

	public  class MyObj
	{
	  public:
			vec_cv3 m_v;
			std::vector < int > m_f;
	  public:
			int get_face(int ifaceNr); //nullbased
			cv3 get_face_normal(int ifaceNr); //nullbased

	};


	public class CGravitation
	{
	public:
		vec_cv3 p, v, f;
		CGravitation();
		void Go();
	};

  public ref class OpenGL_Base
  {
	public:
		int iRender1or2;
	public:
		CGravitation* m_grav;
		MyData m_data;
		MyObj* m_Obj;

	public:
		List<v3^>^ arr_glob1;

		String^ getString(int i);

		List<v3^>^ arr_glob2;
		List<v3^>^ arr_glob3;

  public:
		int m_iDoMouseMove;
		int m_iDoMouseFly;

		double m_speed, m_turn;



		enumScenenAuswahl m_iShowCubes_or_Gel;
		System::Void ToggleShow(enumScenenAuswahl i);

  public:
		HWND m_hWnd1;
		HWND m_hWnd2;
		HDC  m_HDC1;
		HDC  m_HDC2;
		GLint m_iSetPixelFormat1;
		GLint m_iSetPixelFormat2;

		HGLRC m_hglrc1;
		HGLRC m_hglrc2;

    GLfloat	rtri_angle;				// Angle for the triangle
		GLfloat	rquad_angle;				// Angle for the quad
		GLfloat	berg_angle;



    
		GLdouble m_spinX;
    GLdouble m_spinY;
		GLdouble m_transX;
		GLdouble m_transYYY;
		GLdouble m_transZ;

		GLdouble m_spinX0;
    GLdouble m_spinY0;
		GLdouble m_transX0;
    GLdouble m_transY0;

		int m_iMouseX0;
		int m_iMouseY0;

		v3 m_eye, m_look, m_up;





  public:
		OpenGL_Base(System::Windows::Forms::Panel^  panel1, System::Windows::Forms::Panel^  panel2);
    virtual ~OpenGL_Base(void);


		void MyStartOpenGL1(System::Windows::Forms::Panel^  panel);
		void MyStartOpenGL2(System::Windows::Forms::Panel^  panel);

		GLint MySetPixelFormat1(HDC hdc);
		GLint MySetPixelFormat2(HDC hdc);



		bool InitGL(GLvoid);
		GLvoid ReSizeGLScene(int iViewportNr, GLint iSetPixelformat, GLsizei width, GLsizei height);

		System::Void Render_Sphere(System::Void);
		System::Void Render_Planeten(System::Void);
		System::Void Render_Sphere(cv3 p);

		System::Void Render_Wurmloch(System::Void);

		System::Void Render_MinKugeln(System::Void);

		System::Void Render_Sphere_Material(System::Void);

		System::Void Render_House(System::Void);
		System::Void Render_House1(System::Void);
		System::Void Render_Cube(System::Void);
		System::Void Render_Tree(System::Void);
		System::Void Render_Cones(System::Void);
		System::Void Render_WaveBall(System::Void);
		System::Void Render_Gelände(System::Void);
		System::Void Render_Material(System::Void);
		System::Void Render_Material_A(System::Void);
		System::Void Render_Archimedes(System::Void);


		System::Void Render(System::Void);

		System::Void SwapOpenGLBuffers1(System::Void);
		System::Void SwapOpenGLBuffers2(System::Void);

    //KugelTrianguliert

    // Basis-Koerper
#define dd double
    System::Void Basic_Box(dd x, dd y, dd z, dd w, dd l, dd b, dd h, dd& x1, dd& y1);
    System::Void Basic_Disk(dd x, dd y, dd r);
    System::Void Basic_Box_Fenster(dd x, dd y, dd z, dd angle, dd l, dd b, dd h, dd h1, dd h2, dd& x1, dd& y1); //Zwei Boxen übereinander
		void color1();
		void color2();
///////////////////////////
		System::Void ReadObjFile(String^ s);
///////////////////////////

  };



}

