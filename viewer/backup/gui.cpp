///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Mar  6 2014)
// http://www.wxformbuilder.org/
//
// PLEASE DO "NOT" EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#include "gui.h"

///////////////////////////////////////////////////////////////////////////

MyFrame1::MyFrame1( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxFrame( parent, id, title, pos, size, style )
{
	this->SetSizeHints( wxSize( 1100,700 ), wxDefaultSize );
	
	wxBoxSizer* bSizer4;
	bSizer4 = new wxBoxSizer( wxHORIZONTAL );
	
	wxBoxSizer* bSizer1;
	bSizer1 = new wxBoxSizer( wxVERTICAL );
	
	wxBoxSizer* bSizer2;
	bSizer2 = new wxBoxSizer( wxHORIZONTAL );
	
	m_panel1 = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxStaticBoxSizer* sbSizer1;
	sbSizer1 = new wxStaticBoxSizer( new wxStaticBox( m_panel1, wxID_ANY, wxT("harms") ), wxVERTICAL );
	
	m_plot_harms = new MyPlotwindow( this, m_panel1, -1, wxDefaultPosition, wxDefaultSize, wxWANTS_CHARS, wxPLPLOT_BACKEND_GC | wxPLPLOT_DRAW_TEXT );
	m_plot_harms->SetMinSize( wxSize( 400,300 ) );
	
	sbSizer1->Add( m_plot_harms, 1, wxALL|wxEXPAND, 5 );
	
	
	m_panel1->SetSizer( sbSizer1 );
	m_panel1->Layout();
	sbSizer1->Fit( m_panel1 );
	bSizer2->Add( m_panel1, 1, wxALL|wxEXPAND, 5 );
	
	m_panel2 = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxStaticBoxSizer* sbSizer2;
	sbSizer2 = new wxStaticBoxSizer( new wxStaticBox( m_panel2, wxID_ANY, wxT("b-beating") ), wxVERTICAL );
	
	m_plot_bbeats = new MyPlotwindow( this, m_panel2, -1, wxDefaultPosition, wxDefaultSize, wxWANTS_CHARS, wxPLPLOT_BACKEND_GC | wxPLPLOT_DRAW_TEXT );
	m_plot_bbeats->SetMinSize( wxSize( 400,300 ) );
	
	sbSizer2->Add( m_plot_bbeats, 1, wxALL|wxEXPAND, 5 );
	
	
	m_panel2->SetSizer( sbSizer2 );
	m_panel2->Layout();
	sbSizer2->Fit( m_panel2 );
	bSizer2->Add( m_panel2, 1, wxALL|wxEXPAND, 5 );
	
	
	bSizer1->Add( bSizer2, 1, wxEXPAND, 5 );
	
	wxBoxSizer* bSizer3;
	bSizer3 = new wxBoxSizer( wxHORIZONTAL );
	
	m_panel3 = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxStaticBoxSizer* sbSizer3;
	sbSizer3 = new wxStaticBoxSizer( new wxStaticBox( m_panel3, wxID_ANY, wxT("dquads") ), wxVERTICAL );
	
	m_plot_dquads = new MyPlotwindow( this, m_panel3, -1, wxDefaultPosition, wxDefaultSize, wxWANTS_CHARS, wxPLPLOT_BACKEND_GC | wxPLPLOT_DRAW_TEXT );
	m_plot_dquads->SetMinSize( wxSize( 500,300 ) );
	
	sbSizer3->Add( m_plot_dquads, 1, wxALL|wxEXPAND, 5 );
	
	
	m_panel3->SetSizer( sbSizer3 );
	m_panel3->Layout();
	sbSizer3->Fit( m_panel3 );
	bSizer3->Add( m_panel3, 1, wxALL|wxEXPAND, 5 );
	
	m_panel4 = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxStaticBoxSizer* sbSizer4;
	sbSizer4 = new wxStaticBoxSizer( new wxStaticBox( m_panel4, wxID_ANY, wxT("file selection") ), wxVERTICAL );
	
	m_staticText1 = new wxStaticText( m_panel4, wxID_ANY, wxT("label"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText1->Wrap( -1 );
	sbSizer4->Add( m_staticText1, 0, wxALL, 5 );
	
	m_listBox1 = new wxListBox( m_panel4, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0, NULL, 0 ); 
	sbSizer4->Add( m_listBox1, 0, wxALL, 5 );
	
	m_spinCtrl1 = new wxSpinCtrl( m_panel4, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 10, 0 );
	sbSizer4->Add( m_spinCtrl1, 0, wxALL, 5 );
	
	m_spinCtrl2 = new wxSpinCtrl( m_panel4, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 10, 0 );
	sbSizer4->Add( m_spinCtrl2, 0, wxALL, 5 );
	
	m_spinCtrl3 = new wxSpinCtrl( m_panel4, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0, 10, 0 );
	sbSizer4->Add( m_spinCtrl3, 0, wxALL, 5 );
	
	m_textCtrl1 = new wxTextCtrl( m_panel4, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxTE_READONLY );
	sbSizer4->Add( m_textCtrl1, 0, wxALL, 5 );
	
	m_button1 = new wxButton( m_panel4, wxID_ANY, wxT("&Quit"), wxDefaultPosition, wxDefaultSize, 0 );
	sbSizer4->Add( m_button1, 0, wxALL, 5 );
	
	
	m_panel4->SetSizer( sbSizer4 );
	m_panel4->Layout();
	sbSizer4->Fit( m_panel4 );
	bSizer3->Add( m_panel4, 0, wxALL, 5 );
	
	
	bSizer1->Add( bSizer3, 1, wxEXPAND, 5 );
	
	
	bSizer4->Add( bSizer1, 1, wxEXPAND, 5 );
	
	wxBoxSizer* bSizer6;
	bSizer6 = new wxBoxSizer( wxVERTICAL );
	
	m_panel5 = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxStaticBoxSizer* sbSizer5;
	sbSizer5 = new wxStaticBoxSizer( new wxStaticBox( m_panel5, wxID_ANY, wxT("trim quads") ), wxVERTICAL );
	
	m_textCtrl2 = new wxTextCtrl( m_panel5, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE|wxTE_READONLY );
	m_textCtrl2->SetMinSize( wxSize( 200,-1 ) );
	
	sbSizer5->Add( m_textCtrl2, 1, wxALL, 5 );
	
	
	m_panel5->SetSizer( sbSizer5 );
	m_panel5->Layout();
	sbSizer5->Fit( m_panel5 );
	bSizer6->Add( m_panel5, 1, wxALL, 5 );
	
	m_panel6 = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxStaticBoxSizer* sbSizer6;
	sbSizer6 = new wxStaticBoxSizer( new wxStaticBox( m_panel6, wxID_ANY, wxT("weights") ), wxVERTICAL );
	
	m_textCtrl3 = new wxTextCtrl( m_panel6, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxTE_MULTILINE|wxTE_READONLY );
	m_textCtrl3->SetMinSize( wxSize( 200,-1 ) );
	
	sbSizer6->Add( m_textCtrl3, 1, wxALL, 5 );
	
	
	m_panel6->SetSizer( sbSizer6 );
	m_panel6->Layout();
	sbSizer6->Fit( m_panel6 );
	bSizer6->Add( m_panel6, 1, wxALL, 5 );
	
	
	bSizer4->Add( bSizer6, 0, wxEXPAND, 5 );
	
	
	this->SetSizer( bSizer4 );
	this->Layout();
	bSizer4->Fit( this );
	
	this->Centre( wxBOTH );
	
	// Connect Events
	m_listBox1->Connect( wxEVT_COMMAND_LISTBOX_SELECTED, wxCommandEventHandler( MyFrame1::OnChooseLabel ), NULL, this );
	m_button1->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MyFrame1::OnQuitClick ), NULL, this );
}

MyFrame1::~MyFrame1()
{
	// Disconnect Events
	m_listBox1->Disconnect( wxEVT_COMMAND_LISTBOX_SELECTED, wxCommandEventHandler( MyFrame1::OnChooseLabel ), NULL, this );
	m_button1->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MyFrame1::OnQuitClick ), NULL, this );
	
}
