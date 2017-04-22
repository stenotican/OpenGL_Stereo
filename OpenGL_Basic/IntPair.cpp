#include "StdAfx.h"

#include "IntPair.h"
#include <vector>
#include <string>
#include <algorithm> 

CIntPair::CIntPair(void)
{
  i1=0;
  i2=0;
}
CIntPair::CIntPair(int _i1, int _i2)
{
  i1=_i1;
  i2=_i2;
}

CIntPair::~CIntPair(void)
{
}

bool CIntPair::operator<(const CIntPair& ip) const
{
  if(i1<ip.i1)
    return true;
  if(i1>ip.i1)
    return false;
  if(i2<ip.i2)
    return true;
  if(i2>ip.i2)
    return false;
  return false;
}

bool CIntPair::operator==(const CIntPair& ip) const
{
  if(i1 != ip.i1)
    return false;
  if(i2 != ip.i2)
    return false;
  return true;
}

//-------------------------------------------------------------
//-------------------------------------------------------------
//-------------------------------------------------------------

void CIntMatrix::setat(int i1,int i2, int iValue)
{
  CIntPair ip(i1,i2);
  (*this)[ip] = iValue;
}

bool CIntMatrix::find(CIntPair ip) const
{
	bool ret = false;
	std::map<CIntPair, int>::const_iterator it =  std::map<CIntPair, int>::find(ip);
	if (it != end())
	{
		ret = true;
	}
  return ret;
}
bool CIntMatrix::find(int i1, int i2) const
{
  CIntPair ip(i1,i2);
  return find(ip);
}

bool CIntMatrix::find(CIntPair ip, int& iValue) const
{
	bool ret = false;
	std::map<CIntPair, int>::const_iterator it = std::map<CIntPair, int>::find(ip);
	if (it != end())
	{
		iValue = it->second;
		ret = true;
	}
	return ret;
}
//
//
System::String^ CIntMessageMatrix::MessageText(int iMessage)
{
	System::String^ ret = "Unbekannte Message";
  switch(iMessage)
  {
    case 0x0000: ret="WM_NULL"; break;
    case 0x0001: ret="WM_CREATE"; break;
    case 0x0002: ret="WM_DESTROY"; break;
    case 0x0003: ret="WM_MOVE"; break;
    case 0x0005: ret="WM_SIZE"; break;
    case 0x0006: ret="WM_ACTIVATE"; break;
    case 0x0007: ret="WM_SETFOCUS"; break;
    case 0x0008: ret="WM_KILLFOCUS"; break;
    case 0x000A: ret="WM_ENABLE"; break;
    case 0x000B: ret="WM_SETREDRAW"; break;
    case 0x000C: ret="WM_SETTEXT"; break;
    case 0x000D: ret="WM_GETTEXT"; break;
    case 0x000E: ret="WM_GETTEXTLENGTH"; break;
    case 0x000F: ret="WM_PAINT"; break;
    case 0x0010: ret="WM_CLOSE"; break;
    case 0x0011: ret="WM_QUERYENDSESSION"; break;
    case 0x0013: ret="WM_QUERYOPEN"; break;
    case 0x0016: ret="WM_ENDSESSION"; break;
    case 0x0012: ret="WM_QUIT"; break;
    case 0x0014: ret="WM_ERASEBKGND"; break;
    case 0x0015: ret="WM_SYSCOLORCHANGE"; break;
    case 0x0018: ret="WM_SHOWWINDOW"; break;
    case 0x0019: ret="WM_CTLCOLOR"; break;
    case 0x001A: ret="WM_WININICHANGE"; break;
    case 0x001B: ret="WM_DEVMODECHANGE"; break;
    case 0x001C: ret="WM_ACTIVATEAPP"; break;
    case 0x001D: ret="WM_FONTCHANGE"; break;
    case 0x001E: ret="WM_TIMECHANGE"; break;
    case 0x001F: ret="WM_CANCELMODE"; break;
    case 0x0020: ret="WM_SETCURSOR"; break;
    case 0x0021: ret="WM_MOUSEACTIVATE"; break;
    case 0x0022: ret="WM_CHILDACTIVATE"; break;
    case 0x0023: ret="WM_QUEUESYNC"; break;
    case 0x0024: ret="WM_GETMINMAXINFO"; break;
    case 0x0026: ret="WM_PAINTICON"; break;
    case 0x0027: ret="WM_ICONERASEBKGND"; break;
    case 0x0028: ret="WM_NEXTDLGCTL"; break;
    case 0x002A: ret="WM_SPOOLERSTATUS"; break;
    case 0x002B: ret="WM_DRAWITEM"; break;
    case 0x002C: ret="WM_MEASUREITEM"; break;
    case 0x002D: ret="WM_DELETEITEM"; break;
    case 0x002E: ret="WM_VKEYTOITEM"; break;
    case 0x002F: ret="WM_CHARTOITEM"; break;
    case 0x0030: ret="WM_SETFONT"; break;
    case 0x0031: ret="WM_GETFONT"; break;
    case 0x0032: ret="WM_SETHOTKEY"; break;
    case 0x0033: ret="WM_GETHOTKEY"; break;
    case 0x0037: ret="WM_QUERYDRAGICON"; break;
    case 0x0039: ret="WM_COMPAREITEM"; break;
    case 0x003D: ret="WM_GETOBJECT"; break;
    case 0x0041: ret="WM_COMPACTING"; break;
    case 0x0044: ret="WM_COMMNOTIFY"; break;
    case 0x0046: ret="WM_WINDOWPOSCHANGING"; break;
    case 0x0047: ret="WM_WINDOWPOSCHANGED"; break;
    case 0x0048: ret="WM_POWER"; break;
    case 0x004A: ret="WM_COPYDATA"; break;
    case 0x004B: ret="WM_CANCELJOURNAL"; break;
    case 0x004E: ret="WM_NOTIFY"; break;
    case 0x0050: ret="WM_INPUTLANGCHANGEREQUEST"; break;
    case 0x0051: ret="WM_INPUTLANGCHANGE"; break;
    case 0x0052: ret="WM_TCARD"; break;
    case 0x0053: ret="WM_HELP"; break;
    case 0x0054: ret="WM_USERCHANGED"; break;
    case 0x0055: ret="WM_NOTIFYFORMAT"; break;
    case 0x007B: ret="WM_CONTEXTMENU"; break;
    case 0x007C: ret="WM_STYLECHANGING"; break;
    case 0x007D: ret="WM_STYLECHANGED"; break;
    case 0x007E: ret="WM_DISPLAYCHANGE"; break;
    case 0x007F: ret="WM_GETICON"; break;
    case 0x0080: ret="WM_SETICON"; break;
    case 0x0081: ret="WM_NCCREATE"; break;
    case 0x0082: ret="WM_NCDESTROY"; break;
    case 0x0083: ret="WM_NCCALCSIZE"; break;
    case 0x0084: ret="WM_NCHITTEST"; break;
    case 0x0085: ret="WM_NCPAINT"; break;
    case 0x0086: ret="WM_NCACTIVATE"; break;
    case 0x0087: ret="WM_GETDLGCODE"; break;
    case 0x0088: ret="WM_SYNCPAINT"; break;
    case 0x00A0: ret="WM_NCMOUSEMOVE"; break;
    case 0x00A1: ret="WM_NCLBUTTONDOWN"; break;
    case 0x00A2: ret="WM_NCLBUTTONUP"; break;
    case 0x00A3: ret="WM_NCLBUTTONDBLCLK"; break;
    case 0x00A4: ret="WM_NCRBUTTONDOWN"; break;
    case 0x00A5: ret="WM_NCRBUTTONUP"; break;
    case 0x00A6: ret="WM_NCRBUTTONDBLCLK"; break;
    case 0x00A7: ret="WM_NCMBUTTONDOWN"; break;
    case 0x00A8: ret="WM_NCMBUTTONUP"; break;
    case 0x00A9: ret="WM_NCMBUTTONDBLCLK"; break;
    case 0x00AB: ret="WM_NCXBUTTONDOWN"; break;
    case 0x00AC: ret="WM_NCXBUTTONUP"; break;
    case 0x00AD: ret="WM_NCXBUTTONDBLCLK"; break;
    case 0x00FE: ret="WM_INPUT_DEVICE_CHANGE"; break;
    case 0x00FF: ret="WM_INPUT"; break;
    case 0x0100: ret="WM_KEYDOWN"; break;
    case 0x0101: ret="WM_KEYUP"; break;
    case 0x0102: ret="WM_CHAR"; break;
    case 0x0103: ret="WM_DEADCHAR"; break;
    case 0x0104: ret="WM_SYSKEYDOWN"; break;
    case 0x0105: ret="WM_SYSKEYUP"; break;
    case 0x0106: ret="WM_SYSCHAR"; break;
    case 0x0107: ret="WM_SYSDEADCHAR"; break;
    case 0x0109: ret="WM_KEYLAST"; break;
    case 0x0108: ret="WM_KEYLAST"; break;
    case 0x010D: ret="WM_IME_STARTCOMPOSITION"; break;
    case 0x010E: ret="WM_IME_ENDCOMPOSITION"; break;
    case 0x010F: ret="WM_IME_KEYLAST"; break;
    case 0x0110: ret="WM_INITDIALOG"; break;
    case 0x0111: ret="WM_COMMAND"; break;
    case 0x0112: ret="WM_SYSCOMMAND"; break;
    case 0x0113: ret="WM_TIMER"; break;
    case 0x0114: ret="WM_HSCROLL"; break;
    case 0x0115: ret="WM_VSCROLL"; break;
    case 0x0116: ret="WM_INITMENU"; break;
    case 0x0117: ret="WM_INITMENUPOPUP"; break;
    case 0x0119: ret="WM_GESTURE"; break;
    case 0x011A: ret="WM_GESTURENOTIFY"; break;
    case 0x011F: ret="WM_MENUSELECT"; break;
    case 0x0120: ret="WM_MENUCHAR"; break;
    case 0x0121: ret="WM_ENTERIDLE"; break;
    case 0x0122: ret="WM_MENURBUTTONUP"; break;
    case 0x0123: ret="WM_MENUDRAG"; break;
    case 0x0124: ret="WM_MENUGETOBJECT"; break;
    case 0x0125: ret="WM_UNINITMENUPOPUP"; break;
    case 0x0126: ret="WM_MENUCOMMAND"; break;
    case 0x0127: ret="WM_CHANGEUISTATE"; break;
    case 0x0128: ret="WM_UPDATEUISTATE"; break;
    case 0x0129: ret="WM_QUERYUISTATE"; break;
    case 0x0132: ret="WM_CTLCOLORMSGBOX"; break;
    case 0x0133: ret="WM_CTLCOLOREDIT"; break;
    case 0x0134: ret="WM_CTLCOLORLISTBOX"; break;
    case 0x0135: ret="WM_CTLCOLORBTN"; break;
    case 0x0136: ret="WM_CTLCOLORDLG"; break;
    case 0x0137: ret="WM_CTLCOLORSCROLLBAR"; break;
    case 0x0138: ret="WM_CTLCOLORSTATIC"; break;
    case 0x01E1: ret="MN_GETHMENU"; break;
    case 0x0200: ret="WM_MOUSEMOVE"; break;
    case 0x0201: ret="WM_LBUTTONDOWN"; break;
    case 0x0202: ret="WM_LBUTTONUP"; break;
    case 0x0203: ret="WM_LBUTTONDBLCLK"; break;
    case 0x0204: ret="WM_RBUTTONDOWN"; break;
    case 0x0205: ret="WM_RBUTTONUP"; break;
    case 0x0206: ret="WM_RBUTTONDBLCLK"; break;
    case 0x0207: ret="WM_MBUTTONDOWN"; break;
    case 0x0208: ret="WM_MBUTTONUP"; break;
    case 0x0209: ret="WM_MBUTTONDBLCLK"; break;
    case 0x020A: ret="WM_MOUSEWHEEL"; break;
    case 0x020B: ret="WM_XBUTTONDOWN"; break;
    case 0x020C: ret="WM_XBUTTONUP"; break;
    case 0x020D: ret="WM_XBUTTONDBLCLK"; break;
    case 0x0210: ret="WM_PARENTNOTIFY"; break;
    case 0x0211: ret="WM_ENTERMENULOOP"; break;
    case 0x0212: ret="WM_EXITMENULOOP"; break;
    case 0x0213: ret="WM_NEXTMENU"; break;
    case 0x0214: ret="WM_SIZING"; break;
    case 0x0215: ret="WM_CAPTURECHANGED"; break;
    case 0x0216: ret="WM_MOVING"; break;
    case 0x0220: ret="WM_MDICREATE"; break;
    case 0x0221: ret="WM_MDIDESTROY"; break;
    case 0x0222: ret="WM_MDIACTIVATE"; break;
    case 0x0223: ret="WM_MDIRESTORE"; break;
    case 0x0224: ret="WM_MDINEXT"; break;
    case 0x0225: ret="WM_MDIMAXIMIZE"; break;
    case 0x0226: ret="WM_MDITILE"; break;
    case 0x0227: ret="WM_MDICASCADE"; break;
    case 0x0228: ret="WM_MDIICONARRANGE"; break;
    case 0x0229: ret="WM_MDIGETACTIVE"; break;
    case 0x0230: ret="WM_MDISETMENU"; break;
    case 0x0231: ret="WM_ENTERSIZEMOVE"; break;
    case 0x0232: ret="WM_EXITSIZEMOVE"; break;
    case 0x0233: ret="WM_DROPFILES"; break;
    case 0x0234: ret="WM_MDIREFRESHMENU"; break;
    case 0x0281: ret="WM_IME_SETCONTEXT"; break;
    case 0x0282: ret="WM_IME_NOTIFY"; break;
    case 0x0283: ret="WM_IME_CONTROL"; break;
    case 0x0284: ret="WM_IME_COMPOSITIONFULL"; break;
    case 0x0285: ret="WM_IME_SELECT"; break;
    case 0x0286: ret="WM_IME_CHAR"; break;
    case 0x0288: ret="WM_IME_REQUEST"; break;
    case 0x0290: ret="WM_IME_KEYDOWN"; break;
    case 0x0291: ret="WM_IME_KEYUP"; break;
    case 0x02A1: ret="WM_MOUSEHOVER"; break;
    case 0x02A3: ret="WM_MOUSELEAVE"; break;
    case 0x02A0: ret="WM_NCMOUSEHOVER"; break;
    case 0x02A2: ret="WM_NCMOUSELEAVE"; break;
    case 0x02B1: ret="WM_WTSSESSION_CHANGE"; break;
    case 0x02c0: ret="WM_TABLET_FIRST"; break;
    case 0x02df: ret="WM_TABLET_LAST"; break;
    case 0x0300: ret="WM_CUT"; break;
    case 0x0301: ret="WM_COPY"; break;
    case 0x0302: ret="WM_PASTE"; break;
    case 0x0303: ret="WM_CLEAR"; break;
    case 0x0304: ret="WM_UNDO"; break;
    case 0x0305: ret="WM_RENDERFORMAT"; break;
    case 0x0306: ret="WM_RENDERALLFORMATS"; break;
    case 0x0307: ret="WM_DESTROYCLIPBOARD"; break;
    case 0x0308: ret="WM_DRAWCLIPBOARD"; break;
    case 0x0309: ret="WM_PAINTCLIPBOARD"; break;
    case 0x030A: ret="WM_VSCROLLCLIPBOARD"; break;
    case 0x030B: ret="WM_SIZECLIPBOARD"; break;
    case 0x030C: ret="WM_ASKCBFORMATNAME"; break;
    case 0x030D: ret="WM_CHANGECBCHAIN"; break;
    case 0x030E: ret="WM_HSCROLLCLIPBOARD"; break;
    case 0x030F: ret="WM_QUERYNEWPALETTE"; break;
    case 0x0310: ret="WM_PALETTEISCHANGING"; break;
    case 0x0311: ret="WM_PALETTECHANGED"; break;
    case 0x0312: ret="WM_HOTKEY"; break;
    case 0x0317: ret="WM_PRINT"; break;
    case 0x0318: ret="WM_PRINTCLIENT"; break;
    case 0x0319: ret="WM_APPCOMMAND"; break;
    case 0x031A: ret="WM_THEMECHANGED"; break;
    case 0x031D: ret="WM_CLIPBOARDUPDATE"; break;
    case 0x031E: ret="WM_DWMCOMPOSITIONCHANGED"; break;
    case 0x031F: ret="WM_DWMNCRENDERINGCHANGED"; break;
    case 0x0320: ret="WM_DWMCOLORIZATIONCOLORCHANGED"; break;
    case 0x0321: ret="WM_DWMWINDOWMAXIMIZEDCHANGE"; break;
    case 0x0323: ret="WM_DWMSENDICONICTHUMBNAIL"; break;
    case 0x0326: ret="WM_DWMSENDICONICLIVEPREVIEWBITMAP"; break;
    case 0x033F: ret="WM_GETTITLEBARINFOEX"; break;
    case 0x0358: ret="WM_HANDHELDFIRST"; break;
    case 0x035F: ret="WM_HANDHELDLAST"; break;
    case 0x0360: ret="WM_AFXFIRST 360 -----------------"; break;
    case 0x0361: ret="WM_AFXFIRST 361 -----------------"; break;
    case 0x0362: ret="WM_AFXFIRST 362 -----------------"; break;
    case 0x0363: ret="WM_AFXFIRST 363 -----------------"; break;
    case 0x0364: ret="WM_AFXFIRST 364 -----------------"; break;
    case 0x0365: ret="WM_AFXFIRST 365 -----------------"; break;
    case 0x0366: ret="WM_AFXFIRST 366 -----------------"; break;
    case 0x0367: ret="WM_AFXFIRST 367 -----------------"; break;
    case 0x0368: ret="WM_AFXFIRST 368 -----------------"; break;
    case 0x0369: ret="WM_AFXFIRST 369 -----------------"; break;
    case 0x036A: ret="WM_AFXFIRST 36A -----------------"; break;
    case 0x036B: ret="WM_AFXFIRST 36B -----------------"; break;
    case 0x036C: ret="WM_AFXFIRST 36C -----------------"; break;
    case 0x036D: ret="WM_AFXFIRST 36D -----------------"; break;
    case 0x036E: ret="WM_AFXFIRST 36E -----------------"; break;
    case 0x036F: ret="WM_AFXFIRST 36F -----------------"; break;

    case 0x037F: ret="WM_AFXLAST"; break;
    case 0x0380: ret="WM_PENWINFIRST"; break;
    case 0x038F: ret="WM_PENWINLAST"; break;
    case 0x0400: ret="WM_USER_400 ###################"; break;
    case 0x0401: ret="WM_USER_401 ###################"; break;
    case 0x0402: ret="WM_USER_402 ###################"; break;
    case 0x0403: ret="WM_USER_403 ###################"; break;
    case 0x0404: ret="WM_USER_404 ###################"; break;
    case 0x0405: ret="WM_USER_405 ###################"; break;
    case 0x0406: ret="WM_USER_406 ###################"; break;
    case 0x0407: ret="WM_USER_407 ###################"; break;
    case 0x0408: ret="WM_USER_408 ###################"; break;
    case 0x0409: ret="WM_USER_409 ###################"; break;
    case 0x040A: ret="WM_USER_40A ###################"; break;
    case 0x040B: ret="WM_USER_40B ###################"; break;
    case 0x040C: ret="WM_USER_40C ###################"; break;
    case 0x040D: ret="WM_USER_40D ###################"; break;
    case 0x040E: ret="WM_USER_40E ###################"; break;
    case 0x040F: ret="WM_USER_40F ###################"; break;
  }
  return ret;
}
//
//
CDoubleVector::CDoubleVector(const CDoubleVector& objectSrc)
{
	(*this) = objectSrc;
}
void CDoubleVector::operator=(const CDoubleVector& objectSrc)
{
	m_dim = objectSrc.m_dim;
	for(int i=0;i< m_dim;i++)
	{
		this->setat(i,objectSrc.getat(i));
	}
}

//
CDoubleVector::CDoubleVector(int N)
{
	m_dim = N;
	for(int i=0;i<N; i++)
	{
		setat(i,0.0);
	}
}

void CDoubleVector::setat(int ip, double iValue)
{
  (*this)[ip] = iValue;
	//TRACE(_T("%d : %f\n", ip, iValue));
	//System::String^ s;
	//s = s->Format("{0}:{1}", ip, iValue);
	//System::Diagnostics::Trace::WriteLine(s);
}

bool CDoubleVector::find(int ip) const
{
	bool ret = false;
	std::map<int, double>::const_iterator it = std::map<int, double>::find(ip);
	if (it != end())
	{
		//iValue = it->second;
		ret = true;
	}
	return ret;
}


bool CDoubleVector::find(int ip, double& dValue) const
{
	bool ret = false;
	std::map<int, double>::const_iterator it = std::map<int, double>::find(ip);
	if (it != end())
	{
		dValue = it->second;
		ret = true;
	}
	return ret;
}

double CDoubleVector::getat(int i) const
{
	return (*this)(i);
}
void CDoubleVector::add(int i, double d_add)
{
	double d = getat(i);
	d += d_add;
	setat(i,d);
}

void CDoubleVector::add(int NDim, const CDoubleVector& v)
{
  for(int i=0;i<NDim; i++)
	{
		add(i, v.getat(i)); 
	}
}

double CDoubleVector::operator()(int i) const
{
	double ret = 0.0;
	double d=-1;
	if(find(i,d))
		ret = d;
	return ret;
}

void CDoubleVector::swap(int i, int j)
{
	double di = getat(i);
	double dj = getat(j);
	setat(i, dj);
	setat(j, di);
}

double CDoubleVector::sum(int NDim)
{
	double ret = 0.0;
	for (int i = 0; i<NDim; i++)
	{
		ret += getat(i);
	}
	return ret;
}



//
//
//

CDoubleMatrix::CDoubleMatrix(int rows, int cols)
{
	m_nrows = rows;
	m_ncols = cols;
}


void CDoubleMatrix::set_row(int ir, CDoubleVector vec, int Ndim)
{
	for(int ic =0;ic<Ndim; ic++)
	{
		setat(ir,ic,vec(ic));
	}
}

void CDoubleMatrix::get_row(int ir, CDoubleVector& vec, int Ndim)
{
	for(int ic =0;ic<Ndim; ic++)
	{
		vec.setat(ic,getat(ir,ic));
	}
}


double CDoubleMatrix::operator()(int i1,int i2) const
{
	return getat(i1,i2);
}

void CDoubleMatrix::setat(int i1,int i2, double iValue)
{
  CIntPair ip(i1,i2);
  (*this)[ip] = iValue;
}

bool CDoubleMatrix::find(CIntPair ip) const
{
	bool ret = false;
	std::map<CIntPair, double>::const_iterator it = std::map<CIntPair, double>::find(ip);
	if (it != end())
	{
		//iValue = it->second;
		ret = true;
	}
	return ret;
}
bool CDoubleMatrix::find(int i1, int i2) const
{
  CIntPair ip(i1,i2);
  return find(ip);
}

bool CDoubleMatrix::find(CIntPair ip, double& dValue) const
{
	bool ret = false;
	std::map<CIntPair, double>::const_iterator it = std::map<CIntPair, double>::find(ip);
	if (it != end())
	{
		dValue = it->second;
		ret = true;
	}
	return ret;
}

bool CDoubleMatrix::find(int i1, int i2, double& dValue) const
{
  CIntPair ip(i1,i2);
  return find(ip, dValue);
}
double CDoubleMatrix::getat(int i1,int i2) const
{
	double ret = 0.0;
	double d;
	if(find(i1,i2,d))
		ret = d;
	return ret;

}
void CDoubleMatrix::swap(int i1, int i2, int j1, int j2)
{
	double di = getat(i1,i2);
	double dj = getat(j1,j2);
	setat(i1,i2,dj);
	setat(j1,j2,di);
}

//
//
std::vector<std::string> SSTL_Permutation(const char *input) 
{ 
  using namespace std; 

  vector<string> result; 
  string cs = input; 

    stable_sort(cs.begin(), cs.end());   

    do   
        result.insert(result.end(), cs); 
    while(next_permutation(cs.begin(), cs.end())); 

    return result; 
}
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

void test()
{
CDoubleVector v1(4);
CDoubleMatrix m(3,4);
m.setat(1,1,0.4);
}
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
