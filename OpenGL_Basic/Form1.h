#pragma once

#include "OpenGL_Base.h"
#include "SettingsForm.h"

namespace OpenGL_Basic {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	//using namespace System::Xml;
	//using namespace System::Xml::Serialization;
	using namespace System::IO;

	/// <summary>
	/// Summary for Form1
	/// </summary>
	public ref class Form1 : public System::Windows::Forms::Form
	{
	protected:
		/// <summary>
		OpenGLForm::OpenGL_Base ^ m_OpenGL_Base;
		/// </summary>
	public:
		Form1(void)
		{
			InitializeComponent();
			//
			m_OpenGL_Base = gcnew OpenGLForm::OpenGL_Base(panel_main, panel_main2);
			//
		}

	protected:
		/// <summary>
		/// </summary>
		~Form1()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Panel^  panel_main;
	private: System::Windows::Forms::Timer^  timer1;
	private: System::Windows::Forms::MenuStrip^  menuStrip1;
	private: System::Windows::Forms::ToolStripMenuItem^  renderToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  cubesToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  geländeToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  mouseToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  rotateToolStripMenuItem;

	private: System::Windows::Forms::ToolStripMenuItem^  materialToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  sphereMaterialToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  sphereToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  wurmlochToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  kugelnToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  renderToolStripMenuItem1;
	private: System::Windows::Forms::ToolStripMenuItem^  setupToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  dataToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  loadToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  saveToolStripMenuItem;
	private: System::Windows::Forms::OpenFileDialog^  openFileDialog1;
	private: System::Windows::Forms::SaveFileDialog^  saveFileDialog1;

	private: System::Windows::Forms::ToolStripMenuItem^  settingsToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  flyToolStripMenuItem;
	public: System::Windows::Forms::Panel^  panel_main2;
	private: System::Windows::Forms::ToolStripMenuItem^  archimedesToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  fileToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  openToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  saveToolStripMenuItem1;
	private: System::Windows::Forms::ToolStripMenuItem^  planetenToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  cylinderToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  treeToolStripMenuItem;
	public:
	private:


	private: System::ComponentModel::IContainer^  components;
	protected:

	protected:

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>


#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			this->panel_main = (gcnew System::Windows::Forms::Panel());
			this->timer1 = (gcnew System::Windows::Forms::Timer(this->components));
			this->menuStrip1 = (gcnew System::Windows::Forms::MenuStrip());
			this->fileToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->openToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->saveToolStripMenuItem1 = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->renderToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->cubesToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->cylinderToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->planetenToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->geländeToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->materialToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->sphereMaterialToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->sphereToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->wurmlochToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->archimedesToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->mouseToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->rotateToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->flyToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->kugelnToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->renderToolStripMenuItem1 = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->setupToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->dataToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->settingsToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->loadToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->saveToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->openFileDialog1 = (gcnew System::Windows::Forms::OpenFileDialog());
			this->saveFileDialog1 = (gcnew System::Windows::Forms::SaveFileDialog());
			this->panel_main2 = (gcnew System::Windows::Forms::Panel());
			this->treeToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->menuStrip1->SuspendLayout();
			this->SuspendLayout();
			// 
			// panel_main
			// 
			this->panel_main->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->panel_main->BackColor = System::Drawing::SystemColors::ControlDark;
			this->panel_main->Location = System::Drawing::Point(21, 60);
			this->panel_main->Name = L"panel_main";
			this->panel_main->Size = System::Drawing::Size(500, 500);
			this->panel_main->TabIndex = 0;
			this->panel_main->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &Form1::panel_main_Paint);
			this->panel_main->MouseDown += gcnew System::Windows::Forms::MouseEventHandler(this, &Form1::panel_main_MouseDown);
			this->panel_main->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &Form1::panel_main_MouseMove);
			this->panel_main->MouseWheel += gcnew System::Windows::Forms::MouseEventHandler(this, &Form1::panel_main_MouseWheel);
			this->panel_main->Resize += gcnew System::EventHandler(this, &Form1::panel_main_Resize);
			// 
			// timer1
			// 
			this->timer1->Enabled = true;
			this->timer1->Interval = 10;
			this->timer1->Tick += gcnew System::EventHandler(this, &Form1::timer1_Tick);
			// 
			// menuStrip1
			// 
			this->menuStrip1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(5) {
				this->fileToolStripMenuItem,
					this->renderToolStripMenuItem, this->mouseToolStripMenuItem, this->kugelnToolStripMenuItem, this->dataToolStripMenuItem
			});
			this->menuStrip1->Location = System::Drawing::Point(0, 0);
			this->menuStrip1->Name = L"menuStrip1";
			this->menuStrip1->Size = System::Drawing::Size(1105, 24);
			this->menuStrip1->TabIndex = 1;
			this->menuStrip1->Text = L"menuStrip1";
			this->menuStrip1->ItemClicked += gcnew System::Windows::Forms::ToolStripItemClickedEventHandler(this, &Form1::menuStrip1_ItemClicked);
			// 
			// fileToolStripMenuItem
			// 
			this->fileToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {
				this->openToolStripMenuItem,
					this->saveToolStripMenuItem1
			});
			this->fileToolStripMenuItem->Name = L"fileToolStripMenuItem";
			this->fileToolStripMenuItem->Size = System::Drawing::Size(37, 20);
			this->fileToolStripMenuItem->Text = L"File";
			// 
			// openToolStripMenuItem
			// 
			this->openToolStripMenuItem->Name = L"openToolStripMenuItem";
			this->openToolStripMenuItem->Size = System::Drawing::Size(103, 22);
			this->openToolStripMenuItem->Text = L"Open";
			this->openToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::openToolStripMenuItem_Click);
			// 
			// saveToolStripMenuItem1
			// 
			this->saveToolStripMenuItem1->Name = L"saveToolStripMenuItem1";
			this->saveToolStripMenuItem1->Size = System::Drawing::Size(103, 22);
			this->saveToolStripMenuItem1->Text = L"Save";
			// 
			// renderToolStripMenuItem
			// 
			this->renderToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(10) {
				this->cubesToolStripMenuItem,
					this->treeToolStripMenuItem, this->cylinderToolStripMenuItem, this->planetenToolStripMenuItem, this->geländeToolStripMenuItem,
					this->materialToolStripMenuItem, this->sphereMaterialToolStripMenuItem, this->sphereToolStripMenuItem, this->wurmlochToolStripMenuItem,
					this->archimedesToolStripMenuItem
			});
			this->renderToolStripMenuItem->Name = L"renderToolStripMenuItem";
			this->renderToolStripMenuItem->Size = System::Drawing::Size(56, 20);
			this->renderToolStripMenuItem->Text = L"Render";
			// 
			// cubesToolStripMenuItem
			// 
			this->cubesToolStripMenuItem->Checked = true;
			this->cubesToolStripMenuItem->CheckState = System::Windows::Forms::CheckState::Checked;
			this->cubesToolStripMenuItem->Name = L"cubesToolStripMenuItem";
			this->cubesToolStripMenuItem->Size = System::Drawing::Size(158, 22);
			this->cubesToolStripMenuItem->Text = L"Cubes";
			this->cubesToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::cubesToolStripMenuItem_Click);
			// 
			// cylinderToolStripMenuItem
			// 
			this->cylinderToolStripMenuItem->Name = L"cylinderToolStripMenuItem";
			this->cylinderToolStripMenuItem->Size = System::Drawing::Size(158, 22);
			this->cylinderToolStripMenuItem->Text = L"Cylinder";
			this->cylinderToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::cylinderToolStripMenuItem_Click);
			// 
			// planetenToolStripMenuItem
			// 
			this->planetenToolStripMenuItem->Name = L"planetenToolStripMenuItem";
			this->planetenToolStripMenuItem->Size = System::Drawing::Size(158, 22);
			this->planetenToolStripMenuItem->Text = L"Planeten";
			this->planetenToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::planetenToolStripMenuItem_Click);
			// 
			// geländeToolStripMenuItem
			// 
			this->geländeToolStripMenuItem->Name = L"geländeToolStripMenuItem";
			this->geländeToolStripMenuItem->Size = System::Drawing::Size(158, 22);
			this->geländeToolStripMenuItem->Text = L"Gelände";
			this->geländeToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::geländeToolStripMenuItem_Click);
			// 
			// materialToolStripMenuItem
			// 
			this->materialToolStripMenuItem->Name = L"materialToolStripMenuItem";
			this->materialToolStripMenuItem->Size = System::Drawing::Size(158, 22);
			this->materialToolStripMenuItem->Text = L"Material";
			this->materialToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::materialToolStripMenuItem_Click);
			// 
			// sphereMaterialToolStripMenuItem
			// 
			this->sphereMaterialToolStripMenuItem->Name = L"sphereMaterialToolStripMenuItem";
			this->sphereMaterialToolStripMenuItem->Size = System::Drawing::Size(158, 22);
			this->sphereMaterialToolStripMenuItem->Text = L"Sphere_Material";
			this->sphereMaterialToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::sphereMaterialToolStripMenuItem_Click);
			// 
			// sphereToolStripMenuItem
			// 
			this->sphereToolStripMenuItem->Name = L"sphereToolStripMenuItem";
			this->sphereToolStripMenuItem->Size = System::Drawing::Size(158, 22);
			this->sphereToolStripMenuItem->Text = L"Sphere";
			this->sphereToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::sphereToolStripMenuItem_Click);
			// 
			// wurmlochToolStripMenuItem
			// 
			this->wurmlochToolStripMenuItem->Name = L"wurmlochToolStripMenuItem";
			this->wurmlochToolStripMenuItem->Size = System::Drawing::Size(158, 22);
			this->wurmlochToolStripMenuItem->Text = L"Wurmloch";
			this->wurmlochToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::wurmlochToolStripMenuItem_Click);
			// 
			// archimedesToolStripMenuItem
			// 
			this->archimedesToolStripMenuItem->Name = L"archimedesToolStripMenuItem";
			this->archimedesToolStripMenuItem->Size = System::Drawing::Size(158, 22);
			this->archimedesToolStripMenuItem->Text = L"Archimedes";
			this->archimedesToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::archimedesToolStripMenuItem_Click);
			// 
			// mouseToolStripMenuItem
			// 
			this->mouseToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {
				this->rotateToolStripMenuItem,
					this->flyToolStripMenuItem
			});
			this->mouseToolStripMenuItem->Name = L"mouseToolStripMenuItem";
			this->mouseToolStripMenuItem->Size = System::Drawing::Size(55, 20);
			this->mouseToolStripMenuItem->Text = L"Mouse";
			// 
			// rotateToolStripMenuItem
			// 
			this->rotateToolStripMenuItem->Name = L"rotateToolStripMenuItem";
			this->rotateToolStripMenuItem->Size = System::Drawing::Size(108, 22);
			this->rotateToolStripMenuItem->Text = L"Rotate";
			this->rotateToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::rotateToolStripMenuItem_Click);
			// 
			// flyToolStripMenuItem
			// 
			this->flyToolStripMenuItem->Name = L"flyToolStripMenuItem";
			this->flyToolStripMenuItem->Size = System::Drawing::Size(108, 22);
			this->flyToolStripMenuItem->Text = L"Fly";
			this->flyToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::flyToolStripMenuItem_Click);
			// 
			// kugelnToolStripMenuItem
			// 
			this->kugelnToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {
				this->renderToolStripMenuItem1,
					this->setupToolStripMenuItem
			});
			this->kugelnToolStripMenuItem->Name = L"kugelnToolStripMenuItem";
			this->kugelnToolStripMenuItem->Size = System::Drawing::Size(56, 20);
			this->kugelnToolStripMenuItem->Text = L"Kugeln";
			// 
			// renderToolStripMenuItem1
			// 
			this->renderToolStripMenuItem1->Name = L"renderToolStripMenuItem1";
			this->renderToolStripMenuItem1->Size = System::Drawing::Size(111, 22);
			this->renderToolStripMenuItem1->Text = L"Render";
			this->renderToolStripMenuItem1->Click += gcnew System::EventHandler(this, &Form1::renderToolStripMenuItem1_Click);
			// 
			// setupToolStripMenuItem
			// 
			this->setupToolStripMenuItem->Name = L"setupToolStripMenuItem";
			this->setupToolStripMenuItem->Size = System::Drawing::Size(111, 22);
			this->setupToolStripMenuItem->Text = L"Setup";
			// 
			// dataToolStripMenuItem
			// 
			this->dataToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(3) {
				this->settingsToolStripMenuItem,
					this->loadToolStripMenuItem, this->saveToolStripMenuItem
			});
			this->dataToolStripMenuItem->Name = L"dataToolStripMenuItem";
			this->dataToolStripMenuItem->Size = System::Drawing::Size(43, 20);
			this->dataToolStripMenuItem->Text = L"Data";
			// 
			// settingsToolStripMenuItem
			// 
			this->settingsToolStripMenuItem->Name = L"settingsToolStripMenuItem";
			this->settingsToolStripMenuItem->Size = System::Drawing::Size(116, 22);
			this->settingsToolStripMenuItem->Text = L"Settings";
			this->settingsToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::settingsToolStripMenuItem_Click);
			// 
			// loadToolStripMenuItem
			// 
			this->loadToolStripMenuItem->Name = L"loadToolStripMenuItem";
			this->loadToolStripMenuItem->Size = System::Drawing::Size(116, 22);
			this->loadToolStripMenuItem->Text = L"Load";
			// 
			// saveToolStripMenuItem
			// 
			this->saveToolStripMenuItem->Name = L"saveToolStripMenuItem";
			this->saveToolStripMenuItem->Size = System::Drawing::Size(116, 22);
			this->saveToolStripMenuItem->Text = L"Save";
			this->saveToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::saveToolStripMenuItem_Click);
			// 
			// openFileDialog1
			// 
			this->openFileDialog1->FileName = L"openFileDialog1";
			// 
			// panel_main2
			// 
			this->panel_main2->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->panel_main2->BackColor = System::Drawing::SystemColors::ControlDark;
			this->panel_main2->Location = System::Drawing::Point(575, 60);
			this->panel_main2->Name = L"panel_main2";
			this->panel_main2->Size = System::Drawing::Size(500, 500);
			this->panel_main2->TabIndex = 1;
			this->panel_main2->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &Form1::panel_main_Paint2);
			// 
			// treeToolStripMenuItem
			// 
			this->treeToolStripMenuItem->Name = L"treeToolStripMenuItem";
			this->treeToolStripMenuItem->Size = System::Drawing::Size(158, 22);
			this->treeToolStripMenuItem->Text = L"Tree";
			this->treeToolStripMenuItem->Click += gcnew System::EventHandler(this, &Form1::treeToolStripMenuItem_Click);
			// 
			// Form1
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1105, 604);
			this->Controls->Add(this->panel_main2);
			this->Controls->Add(this->panel_main);
			this->Controls->Add(this->menuStrip1);
			this->MainMenuStrip = this->menuStrip1;
			this->Name = L"Form1";
			this->Text = L"Form1";
			this->KeyDown += gcnew System::Windows::Forms::KeyEventHandler(this, &Form1::Form1_KeyDown);
			this->Resize += gcnew System::EventHandler(this, &Form1::Form1_Resize);
			this->menuStrip1->ResumeLayout(false);
			this->menuStrip1->PerformLayout();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void panel_main_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e)
	{
		m_OpenGL_Base->MyStartOpenGL1(panel_main);
		m_OpenGL_Base->Render();
		m_OpenGL_Base->SwapOpenGLBuffers1();
	}
	private: System::Void panel_main_Paint2(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e)
	{
		m_OpenGL_Base->MyStartOpenGL2(panel_main2);
		m_OpenGL_Base->Render();
		m_OpenGL_Base->SwapOpenGLBuffers2();
	}
	private: System::Void panel_main_Resize(System::Object^  sender, System::EventArgs^  e)
	{
		m_OpenGL_Base->ReSizeGLScene(1,m_OpenGL_Base->m_iSetPixelFormat1, panel_main->Width, panel_main->Height);
		m_OpenGL_Base->InitGL();
		Invalidate();
	}
	private: System::Void timer1_Tick(System::Object^  sender, System::EventArgs^  e)
	{
		UNREFERENCED_PARAMETER(sender);
		UNREFERENCED_PARAMETER(e);
		if (m_OpenGL_Base != nullptr)
		{
			m_OpenGL_Base->iRender1or2 = 1;
			m_OpenGL_Base->MyStartOpenGL1(panel_main);
			m_OpenGL_Base->Render();
			m_OpenGL_Base->SwapOpenGLBuffers1();

			m_OpenGL_Base->iRender1or2 = 2;
			m_OpenGL_Base->MyStartOpenGL2(panel_main2);
			m_OpenGL_Base->Render();
			m_OpenGL_Base->SwapOpenGLBuffers2();
			//tb_Info->Text = m_OpenGL_Base->m_data.m_dMinDist.ToString();
		}
	}
	private: System::Void panel_main_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
	{
		m_OpenGL_Base->m_iMouseX0 = e->X;
		m_OpenGL_Base->m_iMouseY0 = e->Y;
		m_OpenGL_Base->m_spinX0 = m_OpenGL_Base->m_spinX;
		m_OpenGL_Base->m_spinY0 = m_OpenGL_Base->m_spinY;
		m_OpenGL_Base->m_transX0 = m_OpenGL_Base->m_transX;
		m_OpenGL_Base->m_transY0 = m_OpenGL_Base->m_transYYY;

		System::String^ s = m_OpenGL_Base->m_spinX0.ToString();
		//tb_Info->Text = s;
	}

	private: System::Void panel_main_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
	{
		if (m_OpenGL_Base->m_iDoMouseMove)
		{
			if (e->Button == ::MouseButtons::Right)
			{
				m_OpenGL_Base->m_transX = 0.005333 * (e->X - m_OpenGL_Base->m_iMouseX0) + m_OpenGL_Base->m_transX0;
				m_OpenGL_Base->m_transYYY = 0.005333 * (e->Y - m_OpenGL_Base->m_iMouseY0) + m_OpenGL_Base->m_transY0;
			}
			else if (e->Button == ::MouseButtons::Left)
			{
				m_OpenGL_Base->m_spinX = 0.5333 * (e->X - m_OpenGL_Base->m_iMouseX0) + m_OpenGL_Base->m_spinX0;
				m_OpenGL_Base->m_spinY = 0.5333 * (e->Y - m_OpenGL_Base->m_iMouseY0) + m_OpenGL_Base->m_spinY0;
			}
			else if (e->Button == ::MouseButtons::Middle)
			{
				m_OpenGL_Base->m_transZ = 0.025333 * (e->Y - m_OpenGL_Base->m_iMouseY0) + m_OpenGL_Base->m_transY0;
			}
			m_OpenGL_Base->ReSizeGLScene(1,m_OpenGL_Base->m_iSetPixelFormat1, panel_main->Width, panel_main->Height);
			m_OpenGL_Base->ReSizeGLScene(2,m_OpenGL_Base->m_iSetPixelFormat2, panel_main->Width, panel_main->Height);
//			m_OpenGL_Base->InitGL();
			Invalidate();
		}

		if (m_OpenGL_Base->m_iDoMouseFly)
		{

		}

	}
					 //
	private: System::Void panel_main_MouseWheel(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
	{
		if (m_OpenGL_Base->m_iDoMouseMove)
		{
			int numberOfTextLinesToMove = e->Delta / 120;
			if (numberOfTextLinesToMove > 0)
			{
				m_OpenGL_Base->m_transZ += 0.1;
			}
			else
			{
				m_OpenGL_Base->m_transZ -= 0.1;
			}
			m_OpenGL_Base->ReSizeGLScene(1,m_OpenGL_Base->m_iSetPixelFormat1,panel_main->Width, panel_main->Height);
			m_OpenGL_Base->ReSizeGLScene(2,m_OpenGL_Base->m_iSetPixelFormat2, panel_main->Width, panel_main->Height);
			//m_OpenGL_Base->InitGL();
			Invalidate();
		}
	}
					 //
	private: System::Void rotateToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
	{
		switch (m_OpenGL_Base->m_iDoMouseMove)
		{
		case 0:
			m_OpenGL_Base->m_iDoMouseMove = 1;
			break;
		default:
			m_OpenGL_Base->m_iDoMouseMove = 0;
			break;
		}
		((System::Windows::Forms::ToolStripMenuItem^)(sender))->Checked = (m_OpenGL_Base->m_iDoMouseMove == 1);
	}
//
private: System::Void flyToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
		switch (m_OpenGL_Base->m_iDoMouseFly)
	{
		case 0:
			m_OpenGL_Base->m_iDoMouseFly = 1;
			break;
		default:
			m_OpenGL_Base->m_iDoMouseFly = 0;
			break;
	}
		((System::Windows::Forms::ToolStripMenuItem^)(sender))->Checked = (m_OpenGL_Base->m_iDoMouseFly == 1);
}

//					 //
	private: System::Void SetCheck_Cubes_or_Gel()
	{
		cubesToolStripMenuItem->Checked = (m_OpenGL_Base->m_iShowCubes_or_Gel == OpenGLForm::eCube);
		treeToolStripMenuItem -> Checked = (m_OpenGL_Base->m_iShowCubes_or_Gel == OpenGLForm::eTree);
		planetenToolStripMenuItem->Checked = (m_OpenGL_Base->m_iShowCubes_or_Gel == OpenGLForm::ePlaneten);
		geländeToolStripMenuItem->Checked = (m_OpenGL_Base->m_iShowCubes_or_Gel == OpenGLForm::eGelände);
		materialToolStripMenuItem->Checked = (m_OpenGL_Base->m_iShowCubes_or_Gel == OpenGLForm::eMaterial);
		sphereMaterialToolStripMenuItem->Checked = (m_OpenGL_Base->m_iShowCubes_or_Gel == OpenGLForm::eSphere_Material);
		sphereToolStripMenuItem->Checked = (m_OpenGL_Base->m_iShowCubes_or_Gel == OpenGLForm::eSphere);
		wurmlochToolStripMenuItem->Checked = (m_OpenGL_Base->m_iShowCubes_or_Gel == OpenGLForm::eWurmloch);
		archimedesToolStripMenuItem->Checked = (m_OpenGL_Base->m_iShowCubes_or_Gel == OpenGLForm::eArchimedes);
	}
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	private: System::Void cubesToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
	{
		m_OpenGL_Base->ToggleShow(OpenGLForm::eCube);
		SetCheck_Cubes_or_Gel();
	}
private: System::Void treeToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
	m_OpenGL_Base->ToggleShow(OpenGLForm::eTree);
	SetCheck_Cubes_or_Gel();
}
private: System::Void planetenToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
	m_OpenGL_Base->ToggleShow(OpenGLForm::ePlaneten);
	SetCheck_Cubes_or_Gel();
}

	private: System::Void geländeToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
	{
		m_OpenGL_Base->ToggleShow(OpenGLForm::eGelände);
		SetCheck_Cubes_or_Gel();
	}


	private: System::Void materialToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
	{
		m_OpenGL_Base->ToggleShow(OpenGLForm::eSphere_Material);
		SetCheck_Cubes_or_Gel();
	}

	private: System::Void sphereMaterialToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
	{
		m_OpenGL_Base->ToggleShow(OpenGLForm::eSphere_Material);
		SetCheck_Cubes_or_Gel();
	}



	private: System::Void sphereToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
	{
		m_OpenGL_Base->ToggleShow(OpenGLForm::eSphere);
		SetCheck_Cubes_or_Gel();
	}
	private: System::Void wurmlochToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
	{
		m_OpenGL_Base->ToggleShow(OpenGLForm::eWurmloch);
		SetCheck_Cubes_or_Gel();
	}

	private: System::Void renderToolStripMenuItem1_Click(System::Object^  sender, System::EventArgs^  e)
	{
		m_OpenGL_Base->ToggleShow(OpenGLForm::eMinKugeln);
	}

		private: System::Void archimedesToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
		{
			m_OpenGL_Base->ToggleShow(OpenGLForm::eArchimedes);
			SetCheck_Cubes_or_Gel();
		}
private: System::Void cylinderToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e) 
{
	m_OpenGL_Base->ToggleShow(OpenGLForm::eCones);
	SetCheck_Cubes_or_Gel();
}

				 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	private: System::Void saveToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
	{
		

		saveFileDialog1->Filter =
			"XML-File|*.xml";
		saveFileDialog1->Title = "Save an XML File";
		saveFileDialog1->ShowDialog();
		// If the file name is not an empty string, open it for saving.
		if (saveFileDialog1->FileName != "")
		{
			System::IO::StreamWriter^ fs = System::IO::File::CreateText(saveFileDialog1->FileName);
			for (int i = 0; i < m_OpenGL_Base->arr_glob1->Count; i++)
			{
				fs->WriteLine(m_OpenGL_Base->getString(i));
			}
			fs->Close();
			//XmlSerializer^ xs = gcnew XmlSerializer(typeof(v3));
			//System::IO::File::CreateText(FileStream fs = new FileStream("Data.xml", FileMode.Create);
			{
				//xs.Serialize(fs, customer);
			}
		}
	}

	private: System::Void settingsToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e) 
	{
		SettingsForm^ sf = gcnew SettingsForm(%m_OpenGL_Base->m_data);
		System::Windows::Forms::DialogResult^ ds = sf->ShowDialog();
		if (*ds == System::Windows::Forms::DialogResult::OK)
		{
			m_OpenGL_Base->arr_glob1 = nullptr;
		}
	}
private: System::Void Form1_Resize(System::Object^  sender, System::EventArgs^  e) 
{
	Control^ control = dynamic_cast<Control^>(sender);
	int iH = control->Size.Height;
	int iW = control->Size.Width;

	int iRandLR = 25;
	int iMitte = 5;

	panel_main->Width = (iW - 2 * iRandLR - iMitte) / 2;
	panel_main2->Width = (iW - 2 * iRandLR - iMitte) / 2;

	panel_main->Height = (iH - 60 - 2*iRandLR);
	panel_main2->Height = (iH - 60 - 2*iRandLR);

	panel_main->Location = Point(iRandLR, 60);
	panel_main2->Location = Point(iRandLR + iMitte + (iW - 2 * iRandLR - iMitte) / 2, 60);

}

	private: System::Void Form1_KeyDown(System::Object^  sender, System::Windows::Forms::KeyEventArgs^  e)
	{
		if (m_OpenGL_Base->m_iDoMouseFly == 1)
		{
			switch (e->KeyCode)
			{
			case Keys::Left:
				break;
			}
		}
  }


private: System::Void openToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e) 
{
	openFileDialog1->DefaultExt = "txt";
	openFileDialog1->Title = "Open Text File";
	openFileDialog1->Filter = "TXT files|*.txt";
	if (openFileDialog1->ShowDialog() == System::Windows::Forms::DialogResult::OK)
		{
			String^path = openFileDialog1->FileName;
			m_OpenGL_Base->ReadObjFile(path);
		}
}


private: System::Void menuStrip1_ItemClicked(System::Object^  sender, System::Windows::Forms::ToolStripItemClickedEventArgs^  e) {
}

};

}//namespace
