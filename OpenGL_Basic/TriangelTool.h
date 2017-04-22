#pragma once
#include <cliext/vector>
#include <gcroot.h>

public ref class punkt
{
public:
  double x,y,z;
public:
  punkt(void){x=y=z=0.0;}
  punkt(const punkt% p){x=p.x;y=p.y;z=p.z;}
  punkt(const punkt^ p){x=p->x;y=p->y;z=p->z;}
};

typedef gcroot<punkt^> wrapped_punkt;

public ref class CTriangelTool
{
public:
  cliext::vector<punkt> m_vP;
public:
  CTriangelTool(void);
};

