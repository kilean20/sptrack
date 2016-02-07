// wxWidgets Application Template

#include "wx/wx.h"
#include "mainframe.h"

class MyApp : public wxApp
{
	public:
		virtual bool OnInit();
};

IMPLEMENT_APP( MyApp )

bool MyApp::OnInit()
{
	MainFrame *frame = new MainFrame(NULL);
	frame->Show();
	return true;
}
