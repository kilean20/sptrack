#include "wx/wx.h"
#include "wxPLplotwindow.h"
#include <vector>
#include <fstream>
#include <wx/txtstrm.h>
#include <wx/wfstream.h>
#include <wx/textfile.h>
#include <wx/regex.h>

class MyPlotwindow : public wxPLplotwindow
{
public:
    MyPlotwindow( wxFrame* frame, wxWindow* parent, wxWindowID id = -1, const wxPoint& pos = wxDefaultPosition,
                  const wxSize& size = wxDefaultSize, long style = 0,
                  int pl_style = wxPLPLOT_NONE );

    void OnChar( wxKeyEvent& event );

private:
    wxFrame* mframe;
};

class PlotParam
{
	public:
		PLFLT ymax;
		PlotParam(){ymax=1.0;}
};

class BarChart
{
	public:
		void plfbox_solid( PLFLT, PLFLT, PLFLT );
		void plfbox_hollow( PLFLT, PLFLT, PLFLT );
	protected:
		wxPLplotstream *pls;
};
class BarChart4 : public BarChart
{
	public:
		BarChart4 (MyPlotwindow*, wxString, PlotParam b4sytle=PlotParam());
};
class BarChart2 : public BarChart
{
	public:
		BarChart2 (MyPlotwindow*, wxString, PlotParam b2sytle=PlotParam());
};
class BbeatPlot 
{
	public:
		BbeatPlot (MyPlotwindow*, wxString, wxString, wxString, PlotParam bbsytle=PlotParam());
	private:
		wxPLplotstream *pls;
		PLFLT ymax;
};
