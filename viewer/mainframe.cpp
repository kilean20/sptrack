#include "mainframe.h"

MainFrame::MainFrame( wxWindow* parent )
:
MyFrame1( parent )
{
	label='b';
	path1="data/harm_tag.dat";
	path20="data/FODO0_tag.twiss";
	path21="data/FODO1_tag.twiss";
	path22="data/FODO2_tag.twiss";
	path3="data/dqba_tag.dat";
	path4="";
	path5="";

	wxTextFile file;
	file.Open("run.sh");
	wxArrayString array;
	//wxRegEx reLabel("(^.)(.+) +(.+?) +(.+?) +(.+?) +(.+?) +(.+?)", wxRE_ADVANCED);
	wxRegEx reLabel("(^.)(\\S+) +(\\S+) +(\\S+) +(\\S+) +(\\S+) +(\\S+)", wxRE_ADVANCED);

	if(file.IsOpened())
	{
		for (wxString str = file.GetFirstLine(); !file.Eof(); str = file.GetNextLine() )
		{
			if(file.GetCurrentLine()>=12)
			{
				reLabel.Matches(str);
				if(reLabel.GetMatch(str,1)=="#")
				{
					array.Add(reLabel.GetMatch(str,2));
					array5.Add(reLabel.GetMatch(str,5));
					array6.Add(reLabel.GetMatch(str,6));
					array7.Add(reLabel.GetMatch(str,7));
				}
			}
		}

	}
	else
	{
		//failed to open! we do nothing for now!
	}

	for(int idx=0; idx<array.GetCount(); idx++)
	{
		m_listBox1->Append(array[idx]);
	}
}

void MainFrame::OnYScaleSet1( wxSpinEvent& event )
{
// TODO: Implement OnYScaleSet1
	label=m_listBox1->GetString(m_listBox1->GetSelection());
	path1.Replace("tag",label,false);
	b4style.ymax = m_spinCtrl1->GetValue();
	BarChart4 *p1 = new BarChart4 ( m_plot_harms, path1, b4style);
	delete p1;
}

void MainFrame::OnYScaleSet2( wxSpinEvent& event )
{
// TODO: Implement OnYScaleSet2
	label=m_listBox1->GetString(m_listBox1->GetSelection());
	path20.Replace("tag",label,false);
	path21.Replace("tag",label,false);
	path22.Replace("tag",label,false);
	bbstyle.ymax = m_spinCtrl2->GetValue();
	BbeatPlot *p2 = new BbeatPlot ( m_plot_bbeats, path20, path21, path22, bbstyle );
	delete p2;
}

void MainFrame::OnYScaleSet3( wxSpinEvent& event )
{
// TODO: Implement OnYScaleSet3
	label=m_listBox1->GetString(m_listBox1->GetSelection());
	path3.Replace("tag",label,false);
	b2style.ymax = m_spinCtrl3->GetValue();
	BarChart2 *p3 = new BarChart2 ( m_plot_dquads, path3, b2style);
	delete p3;
}

void MainFrame::OnChooseLabel( wxCommandEvent& event )
{
// TODO: Implement OnChooseLabel
	label=m_listBox1->GetString(m_listBox1->GetSelection());
std::cout<<label<<std::endl;
	path3.Replace("tag",label,false);
	path20.Replace("tag",label,false);
	path21.Replace("tag",label,false);
	path22.Replace("tag",label,false);
	path1.Replace("tag",label,false);
	BarChart4 *p1 = new BarChart4 ( m_plot_harms, path1, b4style);
	BbeatPlot *p2 = new BbeatPlot ( m_plot_bbeats, path20, path21, path22, bbstyle);
	BarChart2 *p3 = new BarChart2 ( m_plot_dquads, path3, b2style);
	//restore
	path1="data/harm_tag.dat";
	path20="data/FODO0_tag.twiss";
	path21="data/FODO1_tag.twiss";
	path22="data/FODO2_tag.twiss";
	path3="data/dqba_tag.dat";
	delete p1, p2, p3;
	path4=array5[m_listBox1->GetSelection()];
	path5=array6[m_listBox1->GetSelection()];
//std::cout<<path4<<std::endl;

//m_textCtrl1->Show(array7[m_listBox1->GetSelection()]);
	m_textCtrl1->Clear();
	*m_textCtrl1<<array7[m_listBox1->GetSelection()];
	m_textCtrl2->LoadFile(path4);
	m_textCtrl3->LoadFile(path5);
}

void MainFrame::OnQuitClick( wxCommandEvent& event )
{
// TODO: Implement OnSwitchClick
	Close( true );
}



BarChart4::BarChart4(MyPlotwindow *m_plot, wxString filename, PlotParam b4style)
{
    pls = m_plot->GetStream();

	double x[10];
	double xa[10],xb[10],za[10],zb[10];
	wxFileInputStream fin( filename );
	wxTextInputStream text( fin );
	for(int i=0;i<10;i++)
	{
		x[i]=i+1.0;
		text>>xb[i] >>zb[i] >>xa[i] >>za[i];
	}

    //pls->env( 0, 11, 0, 0.012, 0, 0 );
    pls->adv( 0 );
	pls->vsta();
	pls->wind( 0, 11, 0, b4style.ymax*0.015 );
	plcol0( 15 );
    pls->box( "bcnt", 1.0, 0, "bcntv", 0, 4 );
    pls->lab( "harmonic number", "", "Stopband harmonics components" );

	plcol0( 0 );

	plcol1( 1 );
	pls->psty( 0 );
	for (int i = 0; i < 10; i++ )
	{
		plfbox_solid( x[i]-0.4, xb[i], 0.2);
	}
	plcol1( 1 );
	pls->psty( 3 );
	for (int i = 0; i < 10; i++ )
	{
		plfbox_solid( x[i]-0.2, xa[i] , 0.2);
		plfbox_hollow( x[i]-0.2, xa[i] , 0.2);
	}
	plcol1( 0 );
	pls->psty( 0 );
	for (int i = 0; i < 10; i++ )
	{
		plfbox_solid( x[i]+0.0, zb[i], 0.2);
	}
	plcol1( 0 );
	pls->psty( 3 ); //pls->lsty( 1 );
	for (int i = 0; i < 10; i++ )
	{
		plfbox_solid( x[i]+0.2, za[i] , 0.2);
		plfbox_hollow( x[i]+0.2, za[i] , 0.2);
	}

	//legend
	pls->vpor(0.15, 0.90, 0.75, 0.85);
	pls->wind( 0, 1, 0, 1 );
	plcol0( 5 );
	pls->lsty( 1 );
    pls->box( "bc", 0, 0, "bc", 0, 0 );
    pls->ptex( 0.1, 0.70, 0.0, 0.0, 0, "red: horizontal" );
    pls->ptex( 0.1, 0.20, 0.0, 0.0, 0, "blue: vertical" );
    pls->ptex( 0.6, 0.70, 0.0, 0.0, 0, "solid: before" );
    pls->ptex( 0.6, 0.20, 0.0, 0.0, 0, "shadow: after" );
/*
	pls->vpor(0.15, 0.40, 0.65, 0.85);
	pls->wind( 0, 1, 0, 1 );
	plcol0( 5 );
    pls->box( "bc", 0, 0, "bc", 0, 0 );
    pls->ptex( 0.4, 0.6, 0.0, 0.0, 0, "legend 1" );
    pls->ptex( 0.4, 0.3, 0.0, 0.0, 0, "legend 2" );

	pls->vpor(0.6, 0.9, 0.5, 0.85);
	pls->wind( 0, 1, -0.003, 0.003 );
	plcol0( 6 );
    //pls->box( "bc", 0, 0, "bc", 0, 0 );
    pls->box( "abcnt", 0, 0, "bcntv", 0, 0 );
    pls->ptex( 0.5, 0.5, 0.0, 0.0, 0, "inset" );
*/

    m_plot->RenewPlot();
}

BarChart2::BarChart2(MyPlotwindow *m_plot, wxString filename, PlotParam b2style)
{
    pls = m_plot->GetStream();

	wxFileInputStream fin( filename );
	wxTextInputStream text( fin );
	double x[36];
	float a[36],b[36];
	for(int i=0;i<36;i++)
	{
		x[i]=i+1.0;
		text>>b[i]>>a[i];
	}

    pls->adv( 0 );
	pls->vsta();
	pls->wind( 0, 37, -b2style.ymax*0.001, b2style.ymax*0.001 );
	plcol0( 15 );
    pls->box( "abcnt", 5.0, 0, "bcntv", 0, 0 );
    pls->lab( "quad index", "#gDk#d1#ul (1/m)", "" );

	plcol0( 0 );

	plcol1( 1 );
	plpsty( 0 );
	for (int i = 0; i < 36; i++ )
	{
		plfbox_solid( x[i]-0.3, b[i], 0.25);
	}
	plcol1( 0 );
	plpsty( 0 ); //pls->lsty( 1 );
	for (int i = 0; i < 36; i++ )
	{
		plfbox_solid( x[i]+0.05, a[i] , 0.25);
	}

	//legend
	pls->vpor(0.15, 0.28, 0.75, 0.85);
	pls->wind( 0, 1, 0, 1 );
	plcol0( 5 );
    pls->box( "bc", 0, 0, "bc", 0, 0 );
    pls->ptex( 0.4, 0.7, 0.0, 0.0, 0, "before" );
    pls->ptex( 0.4, 0.2, 0.0, 0.0, 0, "after " );

	plcol1( 1 );
	PLFLT x1[]={0.1,0.1,0.35,0.35};
	PLFLT y1[]={0.7,0.8,0.8 ,0.7};
    pls->fill( 4, x1, y1 );

	PLFLT x2[]={0.1,0.1,0.35 ,0.35 };
	PLFLT y2[]={0.2,0.3,0.3  ,0.2};
	plcol1( 0 );
    pls->fill( 4, x2, y2 );

    m_plot->RenewPlot();
}

void BarChart::plfbox_solid( PLFLT x0, PLFLT y0, PLFLT width )
{
    PLFLT *x = new PLFLT[4];
    PLFLT *y = new PLFLT[4];

    x[0] = x0;
    y[0] = 0.;
    x[1] = x0;
    y[1] = y0;
    x[2] = x0 + width;
    y[2] = y0;
    x[3] = x0 + width;
    y[3] = 0.;
    pls->fill( 4, x, y );

    delete[] x;
    delete[] y;
}

void BarChart::plfbox_hollow( PLFLT x0, PLFLT y0, PLFLT width )
{
    PLFLT *x = new PLFLT[4];
    PLFLT *y = new PLFLT[4];

    x[0] = x0;
    y[0] = 0.;
    x[1] = x0;
    y[1] = y0;
    x[2] = x0 + width;
    y[2] = y0;
    x[3] = x0 + width;
    y[3] = 0.;
    pls->line( 4, x, y );

    delete[] x;
    delete[] y;
}

BbeatPlot::BbeatPlot (MyPlotwindow* m_plot, wxString path20, wxString path21, wxString path22, PlotParam bbstyle)
{
    pls = m_plot->GetStream();

	std::vector<PLFLT> S,B1,B2,B3,B4,B5,B6;
	std::ifstream fin1;
	std::ifstream fin2;
	std::ifstream fin3;
	fin1.open(path20.c_str());
	fin2.open(path21.c_str());
	fin3.open(path22.c_str());
	fin1.ignore(1000,'\n');
	fin2.ignore(1000,'\n');
	fin3.ignore(1000,'\n');

	PLFLT s,previous_s,dummy;
	previous_s = -1;
	while(fin1.good())
	{
		fin1.ignore(30,'\n');
		fin2.ignore(30,'\n');
		fin3.ignore(30,'\n');
		fin1>>s;
		fin2.ignore(15,'\n');
		fin3.ignore(15,'\n');
		if(s!=previous_s){
			S.push_back(s);
			fin1.ignore(60,'\n');
			fin2.ignore(60,'\n');
			fin3.ignore(60,'\n');
			fin1>>dummy;
			B1.push_back(dummy);
			fin2>>dummy;
			B2.push_back(dummy);
			fin3>>dummy;
			B3.push_back(dummy);
			fin1.ignore(30,'\n');
			fin2.ignore(30,'\n');
			fin3.ignore(30,'\n');
			fin1>>dummy;
			B4.push_back(dummy);
			fin2>>dummy;
			B5.push_back(dummy);
			fin3>>dummy;
			B6.push_back(dummy);
			fin1.ignore(1000,'\n');
			fin2.ignore(1000,'\n');
			fin3.ignore(1000,'\n');
		}else{
			fin1.ignore(1000,'\n');
			fin2.ignore(1000,'\n');
			fin3.ignore(1000,'\n');
		}
		previous_s = s;
	}
	
	fin1.close();
	fin2.close();
	fin3.close();

	for(int i=0;i<S.size();i++){
		B2[i]=(B2[i]-B1[i])/B1[i]*100;
		B3[i]=(B3[i]-B1[i])/B1[i]*100;
		B5[i]=(B5[i]-B4[i])/B4[i]*100;
		B6[i]=(B6[i]-B4[i])/B4[i]*100;
	}

	PLFLT* x = &S[0];
	PLFLT* y1 = &B2[0];
	PLFLT* y2 = &B3[0];
	PLFLT* y3 = &B5[0];
	PLFLT* y4 = &B6[0];

	pls->adv( 0 );
	pls->col0( 15 );
	pls->env( 0, S.back(), -bbstyle.ymax*2, bbstyle.ymax*2, 0, 0 );
	pls->lab( "s", "(%)", "beat beating" );
	pls->box( "abcnt", 0, 0, "bcntv", 0, 0 );

    pls->width( 1 );

	pls->col0( 1 );
	pls->lsty( 1 );
	pls->line( S.size(), x, y1 );
	pls->col0( 1 );
	pls->lsty( 2 );
	pls->line( S.size(), x, y2 );
	pls->col0( 9 );
	pls->lsty( 1 );
	pls->line( S.size(), x, y3 );
	pls->col0( 9 );
	pls->lsty( 2 );
	pls->line( S.size(), x, y4 );

	//legend
	pls->vpor(0.15, 0.90, 0.75, 0.85);
	pls->wind( 0, 1, 0, 1 );
	plcol0( 5 );
	pls->lsty( 1 );
    pls->box( "bc", 0, 0, "bc", 0, 0 );
    pls->ptex( 0.1, 0.70, 0.0, 0.0, 0, "red: horizontal" );
    pls->ptex( 0.1, 0.20, 0.0, 0.0, 0, "blue: vertical" );
    pls->ptex( 0.6, 0.70, 0.0, 0.0, 0, "solid: before" );
    pls->ptex( 0.6, 0.20, 0.0, 0.0, 0, "dash: after" );
    m_plot->RenewPlot();

}

//! Show information if Menu entry About was choosen.
void MainFrame::OnAbout( wxCommandEvent& WXUNUSED( event ) )
{
	wxMessageBox( _T( "This is the About dialog of the wxPLplot demo.\n" ), _T( "About wxPLplot" ),
			wxOK | wxICON_INFORMATION, this );
}
