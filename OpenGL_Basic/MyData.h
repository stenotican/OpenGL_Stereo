#pragma once


ref class MyData
{
public:
	int m_iAnzahlKugeln;
	int m_iMaxIterations;
	int m_iFormfaktor;
	int m_iOptiZiel; //0= min sum 1/r 1=min dist
	double m_dMinDist;
public:
	MyData();
	virtual ~MyData();
};


