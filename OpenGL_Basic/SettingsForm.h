#pragma once

#include "MyData.h"


namespace OpenGL_Basic {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for SettingsForm
	/// </summary>
	public ref class SettingsForm : public System::Windows::Forms::Form
	{
	public:
		MyData^ m_data;
	public:
		SettingsForm(MyData^ data) :m_data(data)
		{
			InitializeComponent();
			//
			textBox1->Text = m_data->m_iAnzahlKugeln.ToString();
			textBox2->Text = m_data->m_iMaxIterations.ToString();
			textBox3->Text = m_data->m_iFormfaktor.ToString();
      textBox4->Text = m_data->m_iOptiZiel.ToString();
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~SettingsForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::TextBox^  textBox1;
	protected:
	private: System::Windows::Forms::TextBox^  textBox2;
	private: System::Windows::Forms::TextBox^  textBox3;
	private: System::Windows::Forms::Label^  label1;
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::Label^  label3;
	private: System::Windows::Forms::Button^  bt_Ok;
  private: System::Windows::Forms::TextBox^  textBox4;
  private: System::Windows::Forms::Label^  label4;
	private: System::Windows::Forms::TextBox^  textBox_Random;

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->textBox2 = (gcnew System::Windows::Forms::TextBox());
			this->textBox3 = (gcnew System::Windows::Forms::TextBox());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->bt_Ok = (gcnew System::Windows::Forms::Button());
			this->textBox4 = (gcnew System::Windows::Forms::TextBox());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->textBox_Random = (gcnew System::Windows::Forms::TextBox());
			this->SuspendLayout();
			// 
			// textBox1
			// 
			this->textBox1->Location = System::Drawing::Point(108, 44);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(99, 20);
			this->textBox1->TabIndex = 0;
			// 
			// textBox2
			// 
			this->textBox2->Location = System::Drawing::Point(108, 70);
			this->textBox2->Name = L"textBox2";
			this->textBox2->Size = System::Drawing::Size(131, 20);
			this->textBox2->TabIndex = 1;
			// 
			// textBox3
			// 
			this->textBox3->Location = System::Drawing::Point(108, 96);
			this->textBox3->Name = L"textBox3";
			this->textBox3->Size = System::Drawing::Size(48, 20);
			this->textBox3->TabIndex = 2;
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(13, 46);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(64, 13);
			this->label1->TabIndex = 3;
			this->label1->Text = L"Anz Kugeln:";
			this->label1->Click += gcnew System::EventHandler(this, &SettingsForm::label1_Click);
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(12, 77);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(46, 13);
			this->label2->TabIndex = 4;
			this->label2->Text = L"Anz Iter:";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Location = System::Drawing::Point(13, 103);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(60, 13);
			this->label3->TabIndex = 5;
			this->label3->Text = L"Formfaktor:";
			// 
			// bt_Ok
			// 
			this->bt_Ok->DialogResult = System::Windows::Forms::DialogResult::OK;
			this->bt_Ok->Location = System::Drawing::Point(362, 134);
			this->bt_Ok->Name = L"bt_Ok";
			this->bt_Ok->Size = System::Drawing::Size(94, 30);
			this->bt_Ok->TabIndex = 6;
			this->bt_Ok->Text = L"Ok";
			this->bt_Ok->UseVisualStyleBackColor = true;
			this->bt_Ok->Click += gcnew System::EventHandler(this, &SettingsForm::bt_Ok_Click);
			// 
			// textBox4
			// 
			this->textBox4->Location = System::Drawing::Point(108, 134);
			this->textBox4->Name = L"textBox4";
			this->textBox4->Size = System::Drawing::Size(48, 20);
			this->textBox4->TabIndex = 7;
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Location = System::Drawing::Point(13, 141);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(42, 13);
			this->label4->TabIndex = 8;
			this->label4->Text = L"OptiArt:";
			// 
			// textBox_Random
			// 
			this->textBox_Random->Location = System::Drawing::Point(191, 96);
			this->textBox_Random->Name = L"textBox_Random";
			this->textBox_Random->Size = System::Drawing::Size(48, 20);
			this->textBox_Random->TabIndex = 9;
			// 
			// SettingsForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(471, 176);
			this->Controls->Add(this->textBox_Random);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->textBox4);
			this->Controls->Add(this->bt_Ok);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->textBox3);
			this->Controls->Add(this->textBox2);
			this->Controls->Add(this->textBox1);
			this->Name = L"SettingsForm";
			this->Text = L"SettingsForm";
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void label1_Click(System::Object^  sender, System::EventArgs^  e) {
	}

	private: System::Void bt_Ok_Click(System::Object^  sender, System::EventArgs^  e) 
	{
		String^s = textBox1->Text;
		int i1,i2;
		i2 = i1.Parse(s);
		
		int i3 = int::Parse(textBox1->Text);
		m_data->m_iAnzahlKugeln = i3;
		m_data->m_iMaxIterations = int::Parse(textBox2->Text);
		m_data->m_iFormfaktor = int::Parse(textBox3->Text);
    m_data->m_iOptiZiel = int::Parse(textBox4->Text);
		Close();
	}
};
}
