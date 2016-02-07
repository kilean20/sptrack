#include "custom.h"
MyPlotwindow::MyPlotwindow( wxFrame* frame, wxWindow* parent, wxWindowID id, const wxPoint& pos,
                            const wxSize& size, long style, int pl_style ) :
    wxPLplotwindow( parent, id, pos, size, style, pl_style )
{
    mframe = frame;
}

void MyPlotwindow::OnChar( wxKeyEvent& event )
{
    int keycode = event.GetKeyCode();

    if ( keycode == WXK_RETURN ||
         keycode == WXK_SPACE ||
         keycode == WXK_RIGHT ||
         keycode == WXK_ESCAPE )
        mframe->Close( true );
    else
        event.Skip();
}
