///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Mar  6 2014)
// http://www.wxformbuilder.org/
//
// PLEASE DO "NOT" EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#ifndef __GUI_H__
#define __GUI_H__

#include <wx/artprov.h>
#include <wx/xrc/xmlres.h>
#include "custom.h"
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/string.h>
#include <wx/stattext.h>
#include <wx/spinctrl.h>
#include <wx/sizer.h>
#include <wx/statbox.h>
#include <wx/panel.h>
#include <wx/listbox.h>
#include <wx/textctrl.h>
#include <wx/button.h>
#include <wx/frame.h>

///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
/// Class MyFrame1
///////////////////////////////////////////////////////////////////////////////
class MyFrame1 : public wxFrame 
{
	private:
	
	protected:
		wxPanel* m_panel1;
		MyPlotwindow* m_plot_harms;
		wxStaticText* m_staticText2;
		wxSpinCtrl* m_spinCtrl1;
		wxPanel* m_panel2;
		MyPlotwindow* m_plot_bbeats;
		wxStaticText* m_staticText3;
		wxSpinCtrl* m_spinCtrl2;
		wxPanel* m_panel3;
		MyPlotwindow* m_plot_dquads;
		wxStaticText* m_staticText4;
		wxSpinCtrl* m_spinCtrl3;
		wxPanel* m_panel4;
		wxStaticText* m_staticText1;
		wxListBox* m_listBox1;
		wxTextCtrl* m_textCtrl1;
		wxButton* m_button1;
		wxPanel* m_panel5;
		wxTextCtrl* m_textCtrl2;
		wxPanel* m_panel6;
		wxTextCtrl* m_textCtrl3;
		
		// Virtual event handlers, overide them in your derived class
		virtual void OnYScaleSet1( wxSpinEvent& event ) { event.Skip(); }
		virtual void OnYScaleSet2( wxSpinEvent& event ) { event.Skip(); }
		virtual void OnYScaleSet3( wxSpinEvent& event ) { event.Skip(); }
		virtual void OnChooseLabel( wxCommandEvent& event ) { event.Skip(); }
		virtual void OnQuitClick( wxCommandEvent& event ) { event.Skip(); }
		
	
	public:
		
		MyFrame1( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("Data Viewer"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( -1,-1 ), long style = wxDEFAULT_FRAME_STYLE|wxTAB_TRAVERSAL );
		
		~MyFrame1();
	
};

#endif //__GUI_H__
