/******************************************************************
 * Modul name : mrf.cpp
 * Author     : Csaba Gradwohl (Gradwohl.Csaba@stud.u-szeged.hu) with some
 *              minor contributions from Zoltan Kato (kato@inf.u-szeged.hu).
 * Copyright  : GNU General Public License www.gnu.org/copyleft/gpl.html
 * Description:
 * Intensity-based image segmentation using a Markov random field
 * segmentation model and four different optimization algorithms:
 * Metropolis - Simulated Annealing using Metropolis dynamics
 * Gibbs      - Simulated Annealing using a Gibbs sampler
 * ICM        - Iterated Conditional Modes, a deterministic suboptimal
 *              method (depends on a good initialization).
 * MMD        - Modified Metropolis Dynamics, a pseudo-stochastic
 *              suboptimal method which is less sensitive to
 *              initialization than ICM.
 *
 * The program GUI is written in wxWindows hence the code can be
 * compiled and ran under Windows as well as under Linux/Unix.
 *
 * $Id: mrf.cpp,v 1.8 2009/01/08 15:30:47 kato Exp $
 * $Revision: 1.8 $
 * $State: Exp $
 * $Log: mrf.cpp,v $
 * Revision 1.8  2009/01/08 15:30:47  kato
 * Upgraded for wxWidgets 2.6.x and later.
 *
 * Revision 1.7  2005/02/15 21:46:16  kato
 * CPU timer now doesn't count GUI overhead
 *
 * Revision 1.6  2005/02/14 15:25:38  kato
 * fixed frame size
 *
 * Revision 1.5  2004/12/13 13:35:04  kato
 * Adjusted position of version & copyright info
 *
 * Revision 1.4  2004/12/13 13:31:12  kato
 * Added copyright & version information,
 * Class parameter window is larger
 *
 * Revision 1.3  2004/12/13 13:15:02  kato
 * Bug fixes
 *
 * Revision 1.2  2004/12/13 12:07:28  kato
 * added CPU timer
 *
 * Revision 1.1  2004/12/08 13:15:06  kato
 * Initial revision
 *
 * 
 *****************************************************************/
#ifndef lint
static char rcsid_mrf_cpp[]="$Id: mrf.cpp,v 1.8 2009/01/08 15:30:47 kato Exp $";
#endif

/* wxWindows includes
 */
#include <wx/wxprec.h>
#ifndef WX_PRECOMP
    #include <wx/wx.h>
#endif
#include <wx/image.h>

#include <math.h>
#include <stdlib.h>

/* Random number generators
 */
#include "randomc.h"   // define classes for random number generators

/* Timer classes
 */
#include "CKProcessTimeCounter.h"

#define WINDOW_TITLE "MRF Image Segmentation Demo $Revision: 1.8 $"
#define VERSION      "MRF Image Segmentation Demo $Revision: 1.8 $ (Last built "\
                     __DATE__" "__TIME__") "
#define COPYRIGHT    "(c) 2004 by Csaba Gradwohl & Zoltan Kato (SZTE - Hungary)"


static wxTextCtrl *gaussians;            // output textfield for Gaussian parameters
static CKProcessTimeCounter timer("core"); // CPU timer
static bool timer_valid = FALSE;

/* Program's application class 
 */
class MyApp: public wxApp
{
  virtual bool OnInit(); // this is the main entry point
};


/* ImageOperations class: it handles all image operations such as
 * loading, saving, etc... 
 */
class ImageOperations
{
public:
  ImageOperations(wxWindow *_frame);    // constructor
  wxImage *LoadBmp(wxString bmp_name);	// loads an image from file		
  bool SaveBmp(wxString bmp_name);      // saves out_image to a given file
  bool IsOutput();			// TRUE if  out_image <> NULL
  void SetNoRegions(int n);	      	// sets the number of regions,
					// allocates/frees memory for
					// means and variances
  int GetNoRegions() { return no_regions; }
  void SetBeta(double b) { beta = b; }
  void SetT(double x) { t = x; }
  void SetT0(double t) { T0 = t; }
  void SetC(double x) { c = x; }
  void SetAlpha(double x) { alpha = x; }
  int GetK() { return K; }
  double GetT() { return T; }
  double GetE() { return E; }
  double GetTimer() { return (timer_valid? timer.GetElapsedTimeMs() : 0.0); }

  void CalculateMeanAndVariance(int region);  // computes mean and
					      // variance of the given region.
  double CalculateEnergy();                   // computes global energy
					      // based on the current
					      // lableing in data
  double LocalEnergy(int i, int j, int label);// computes the local
					      // energy at site (i,j)
					      // assuming "label" has
					      // been assigned to it.

  void Metropolis(bool mmd=false);  // executes Metropolis or MMD (if mmd=true)
  void ICM();			    // executes ICM
  void Gibbs();			    // executes Gibbs sampler

private:
  wxWindow *frame;		    // the main window
  wxImage *in_image, *out_image;    // input & output images
  int width, height;		    // width and height of the
				    // displayed image
  int no_regions;	            // number of regions for Gaussian
				    // parameter computation
  double beta;                      // strength of second order clique potential
  double t;			    // Stop criteraia threshold: stop
				    // if (deltaE < t)
  double T0;		            // Initial temperature (not used by ICM)

  double c;			    // Temperature scheduler's factor:
				    // T(n+1)=c*T(n).
  double alpha;		            // alpha value for MMD
  double *mean;			    // computed mean values and
  double *variance;		    // variances for each region
  double E;			    // current global energy
  double E_old;			    // global energy in the prvious iteration
  double T;			    // current temperature
  int K;			    // current iteration #
  int **classes;		    // this is the labeled image
  int **in_image_data;		    // Intensity values of the input image 

  void InitOutImage();
  void CreateOutput();	           // creates and draws the output
				   // image based on the current labeling
  double Singleton(int i, int j, int label); // computes singleton
					     // potential at site
					     // (i,j) having a label "label"
  double Doubleton(int i, int j, int label); // computes doubleton
					     // potential at site
					     // (i,j) having a label "label"
};


/* MyScrolledWindow class: the window used for diaplaying images
 */
class MyScrolledWindow: public wxScrolledWindow
{
public:
  MyScrolledWindow(wxWindow* parent, wxWindowID id = -1, 
		   const wxPoint& pos = wxDefaultPosition, 
		   const wxSize& size = wxDefaultSize, 
		   long style = wxHSCROLL | wxVSCROLL, 
		   const wxString& name = "scrolledWindow"): 
    wxScrolledWindow(parent, id, pos, size, style, name) 
  { bmp = NULL; }
  void SetBmp(wxImage *_bmp);     // assigns the image to the window.

protected:
  virtual void OnDraw(wxDC& dc);  // displays the image in the window

private:
  wxImage *bmp;			  // the image to be displayed
  int xDst, yDst;		  // the position of the image within
				  // the window (meaningful only when
				  // the image is smaller than the window)

  wxMemoryDC memDC;		  // memDC storing the image

  void OnLeftDown(wxMouseEvent& event);	   // Left button event handler
  void OnMouseMotion(wxMouseEvent& event); // mouse motion event handler
  DECLARE_EVENT_TABLE()
};


/* MyFrame class: the main window
 */
class MyFrame: public wxFrame
{
public:
  MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size);
  ~MyFrame();
  /* returns the coordinates of the given training rectangle (or the
   * current one if region==-1)
   */
  void GetRegion(int &x, int &y, int &w, int &h, int region=-1);
  
  int GetActRegion() { return act_region;}
  void SetRegs1(int x, int y) { 
    regs[act_region*4] = x; 
    regs[act_region*4+1] = y; 
  }
  void SetRegs2(int x, int y) { 
    regs[act_region*4+2] = x; 
    regs[act_region*4+3] = y; 
  }
  MyScrolledWindow *GetInputWindow() { return input_window; }
  MyScrolledWindow *GetOutputWindow() { return output_window; }
  
  bool IsSelected(int region) { // tells whether the region has been selected 
    return (regs[region*4] || regs[region*4+1] || 
	    regs[region*4+2] || regs[region*4+3]); 
  } 
  bool AllRegionsSelected() { // tells whether all the regions has been selected
    for(int i=1; i<imageop->GetNoRegions(); ++i)
      if (!IsSelected(i)) return false; 
    return true; 
  } 
  
private:
  ImageOperations *imageop; 
  MyScrolledWindow *input_window, *output_window;    // input & output
						     // images' window
  wxButton *load_button, *save_button, *doit_button; // buttons
  wxButton *select_region_button;
  wxChoice *op_choice;		// scroll-list of optimization algorithms
  wxTextCtrl *regions;          // input field for number of classes,
  wxTextCtrl *tbeta, *tt;	// beta, threshold t,
  wxTextCtrl *tT0, *tc;		// initial temperature T0, scheduler factor c,
  wxTextCtrl *talpha;		// and MMD's alpha
  int act_region;   // the current class
  int *regs;	    // stores the training rectangles for each class.
  
  /* Event handlers
   */  
  void OnOpen(wxCommandEvent& event);         // Load
  void OnSave(wxCommandEvent& event);         // Save
  void OnDoit(wxCommandEvent& event);         // DoIt
  void OnChoice(wxCommandEvent& event);	      // optimization method selection 
  void OnRegions(wxCommandEvent& event);      // number of classes
  void OnSelectRegion(wxCommandEvent& event); // select training rectangle
  void OnPaint(wxPaintEvent& event);	      // paint handler
  DECLARE_EVENT_TABLE()
    };

enum { ID_LOAD_BUTTON, ID_SAVE_BUTTON, ID_DOIT_BUTTON, ID_CHOICE,
       ID_REGIONS, ID_SELECTREGION_BUTTON, ID_BETA, ID_T, ID_T0, ID_C,
       ID_ALPHA, ID_GAUSSIANS };

/* Event table
 */
BEGIN_EVENT_TABLE(MyFrame, wxFrame)
  EVT_BUTTON(ID_LOAD_BUTTON, MyFrame::OnOpen)
  EVT_BUTTON(ID_SAVE_BUTTON, MyFrame::OnSave)
  EVT_BUTTON(ID_DOIT_BUTTON, MyFrame::OnDoit)
  EVT_CHOICE(ID_CHOICE, MyFrame::OnChoice)
  EVT_PAINT(MyFrame::OnPaint)
  EVT_TEXT(ID_REGIONS, MyFrame::OnRegions)
  EVT_BUTTON(ID_SELECTREGION_BUTTON, MyFrame::OnSelectRegion)
END_EVENT_TABLE()

BEGIN_EVENT_TABLE(MyScrolledWindow, wxScrolledWindow)
  EVT_LEFT_DOWN(MyScrolledWindow::OnLeftDown)
  EVT_MOTION(MyScrolledWindow::OnMouseMotion)
END_EVENT_TABLE()

IMPLEMENT_APP(MyApp)

bool MyApp::OnInit()
{
  MyFrame *frame = new MyFrame( WINDOW_TITLE, 
				wxPoint(0,0), wxSize(900,620) );
  frame->Show( TRUE );
  SetTopWindow( frame );
  return TRUE;
}


void MyScrolledWindow::SetBmp(wxImage *_bmp) 
{
  xDst = yDst = 0;
  // center if image is smaller than the window
  if (_bmp!=NULL)
    {
      if (_bmp->GetWidth() < 300) xDst = (300-_bmp->GetWidth())/2;
      if (_bmp->GetHeight() < 250) yDst = (250-_bmp->GetHeight())/2;
#if wxCHECK_VERSION(2,6,0) // for version 2.6.0 or later
	  wxBitmap *mp=new wxBitmap((const wxImage&)*_bmp);
	  memDC.SelectObject((wxBitmap&)*mp); 
#else    // for version 2.4.x
      memDC.SelectObject(*_bmp);
#endif
    }
  bmp = _bmp; 
}

	
void MyScrolledWindow::OnDraw(wxDC& dc)
{
  if (bmp != NULL)
    {
      // determine which part of the image is visible in the window.
      int x, y;
      GetViewStart(&x, &y);
      x *= 10; y *= 10;	   // must be multiplied by ScrollUnit
      // copy the visible part into the window
      int _xDst, _yDst;
      CalcUnscrolledPosition(xDst, yDst, &_xDst, &_yDst);
      dc.Blit(_xDst, _yDst, 300, 250, &memDC, x, y);
      
      // draw the training rectangle on the input image
      MyFrame *parent = (MyFrame *)GetParent();
      if (parent->GetInputWindow() == this)
	{
	  int x1, y1, w, h;
	  parent->GetRegion(x1, y1, w, h);
	  if (x1!=0 || y1!=0 || w!=0 || h!=0)
	    {
	      wxPen pen(*wxRED_PEN);
	      wxBrush brush(*wxTRANSPARENT_BRUSH);
	      dc.SetPen(pen);
	      dc.SetBrush(brush);
	      dc.DrawRectangle(x1+xDst, y1+yDst, w, h);
	    }
	}
    }
}


/* Left mouse button event handler
 */
void MyScrolledWindow::OnLeftDown(wxMouseEvent& event)
{
  // comvert window coordinates to image coordinates
  int x, y;
  GetViewStart(&x, &y);
  x *= 10; y *= 10;	 // must be multiplied by ScrollUnit

  MyFrame *frame = (MyFrame *)GetParent();
  if (frame->GetActRegion() != -1)	// in this case regs != NULL
    {
      if (event.m_x >= xDst && event.m_x < xDst+bmp->GetWidth() &&
	  event.m_y >= yDst && event.m_y < yDst+bmp->GetHeight())
	frame->SetRegs1(event.m_x+x-xDst, event.m_y+y-yDst); // scroll added
	}
}


/* Mouse motion event handler
 */
void MyScrolledWindow::OnMouseMotion(wxMouseEvent& event)
{
  if (event.LeftIsDown())
    {
      // comvert window coordinates to image coordinates
      int x, y;
      GetViewStart(&x, &y);
      x *= 10; y *= 10;		// must be multiplied by ScrollUnit
      MyFrame *frame = (MyFrame *)GetParent();
      if (frame->GetActRegion() != -1)		// in this case regs != NULL
	{
	  if (event.m_x >= xDst && event.m_x < xDst+bmp->GetWidth() &&
	      event.m_y >= yDst && event.m_y < yDst+bmp->GetHeight())
	    frame->SetRegs2(event.m_x+x-xDst, event.m_y+y-yDst); // scroll added
	  Refresh();
		}
	}
}


/*********************************************************************
/* Functions of MyFrame class
/********************************************************************/
void MyFrame::OnPaint(wxPaintEvent& event)
{
  wxPaintDC pDC(this);

  wxString str;
  str.Printf("Number of classes:");
  pDC.DrawText(str, 20, 325);
  
  str.Printf("ß = ");
  pDC.DrawText(str, 43, 360);
  str.Printf("t = ");
  pDC.DrawText(str, 190, 360);
  str.Printf("Class parameters:");
  pDC.DrawText(str, 20, 465);
  if (op_choice->GetStringSelection() != "ICM") 
    {
      str.Printf("T0 = ");
      pDC.DrawText(str, 35, 395);
      str.Printf("c = ");
      pDC.DrawText(str, 190, 395);
    }
  str.Printf(VERSION);
  pDC.DrawText(str, 330, 550);
  str.Printf(COPYRIGHT);
  pDC.DrawText(str, 330, 565);
  

  if (op_choice->GetStringSelection() == "MMD") 
    {
      str.Printf("alpha = ");
      pDC.DrawText(str, 15, 430);
    }
  
  str.Printf("iteration = ");
  pDC.DrawText(str, 535, 360);
  str.Printf("global energy = ");
  pDC.DrawText(str, 535, 395);
  str.Printf("T = ");
  pDC.DrawText(str, 535, 430);
  str.Printf("CPU time = ");
  pDC.DrawText(str, 535, 465);
  pDC.DrawText(wxString() << imageop->GetK(), 645, 360);
  pDC.DrawText(wxString() << imageop->GetE(), 645, 395);
  pDC.DrawText(wxString() << imageop->GetT(), 645, 430);
  pDC.DrawText(wxString() << imageop->GetTimer() << " ms", 645, 465);
  event.Skip();	
}


MyFrame::MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
  : wxFrame((wxFrame *)NULL, -1, title, pos, size)
{
  imageop = new ImageOperations(this);
  input_window = new MyScrolledWindow(this, -1, wxPoint(15,15), 
				      wxSize(300,250), wxHSCROLL | wxVSCROLL);
  output_window = new MyScrolledWindow(this, -1, wxPoint(475,15), 
				       wxSize(300,250), wxHSCROLL | wxVSCROLL);
  load_button = new wxButton(this, ID_LOAD_BUTTON, "Load", wxPoint(125,280));
  save_button = new wxButton(this, ID_SAVE_BUTTON, "Save", wxPoint(595,280));
  doit_button = new wxButton(this, ID_DOIT_BUTTON, "Do it >>", 
			     wxPoint(358,150));
  doit_button->Disable();
  select_region_button = new wxButton(this, ID_SELECTREGION_BUTTON, 
				      "Select classes", wxPoint(218,321));
  select_region_button->Disable();
  wxString choices[4] = {"Metropolis", "Gibbs sampler", "ICM", "MMD"};
  op_choice = new wxChoice(this, ID_CHOICE, wxPoint(346,80), wxDefaultSize, 
			   4, choices);
  op_choice->SetStringSelection("Metropolis");
	
  regions = new wxTextCtrl(this, ID_REGIONS, "", wxPoint(152,321), 
			   wxSize(27,20), 
			   wxTE_PROCESS_ENTER|wxTE_RIGHT, 
			   *(new wxTextValidator(wxFILTER_NUMERIC)));
  regions->SetMaxLength(3);
  regions->Disable();
  tbeta = new wxTextCtrl(this, ID_BETA, "", wxPoint(67,356), 
			 wxSize(60,20), 
			 wxTE_PROCESS_ENTER|wxTE_RIGHT, 
			 *(new wxTextValidator(wxFILTER_NUMERIC)));
  tbeta->SetMaxLength(8);
  //  tbeta->SetValue("0.9");
  *tbeta << 0.9;
  tt = new wxTextCtrl(this, ID_T, "", wxPoint(212,356), 
		      wxSize(60,20), wxTE_PROCESS_ENTER|wxTE_RIGHT, 
		      *(new wxTextValidator(wxFILTER_NUMERIC)));
  tt->SetMaxLength(8);
  // tt->SetValue("0.05");
  *tt << 0.05;
  tT0 = new wxTextCtrl(this, ID_T0, "", wxPoint(67,391), wxSize(60,20), 
		       wxTE_PROCESS_ENTER|wxTE_RIGHT, 
		       *(new wxTextValidator(wxFILTER_NUMERIC)));
  tT0->SetMaxLength(8);
  // tT0->SetValue("4.0");
  *tT0 << 4.0;
  tc = new wxTextCtrl(this, ID_C, "", wxPoint(212,391), wxSize(60,20), 
		      wxTE_PROCESS_ENTER|wxTE_RIGHT, 
		      *(new wxTextValidator(wxFILTER_NUMERIC)));
  tc->SetMaxLength(8);
  // tc->SetValue("0.98");
  *tc << 0.98;
  talpha = new wxTextCtrl(this, ID_ALPHA, "", wxPoint(67,426), 
			 wxSize(60,20), wxTE_PROCESS_ENTER|wxTE_RIGHT, 
			 *(new wxTextValidator(wxFILTER_NUMERIC)));
  talpha->SetMaxLength(8);
  // talpha->SetValue("0.1");
  *talpha << 0.1;
  talpha->Hide();

  gaussians = new wxTextCtrl(this, ID_GAUSSIANS, "", wxPoint(20,480), 
			     wxSize(300,100), 
			     wxTE_MULTILINE|wxTE_DONTWRAP|wxTE_READONLY,
			     wxDefaultValidator);
  gaussians->SetValue("#\tMean\t\tVariance\n");

  regs = NULL;
  act_region = -1;

}


MyFrame::~MyFrame()
{
  delete imageop;
}


void MyFrame::OnOpen(wxCommandEvent& event)
{
  wxString image_name;
  wxFileDialog* fdialog = new wxFileDialog(this, "Open file", "", "", 
					   "BMP files (*.bmp)|*.bmp", 
					   wxOPEN|wxCHANGE_DIR);
	
  if (fdialog->ShowModal() == wxID_OK)
    {
      image_name = fdialog->GetPath();
      wxImage *bmp;
      if (bmp = imageop->LoadBmp(image_name)) // image succesfully load
	{
	  input_window->SetScrollbars(10,10,(bmp->GetWidth())/10,
				      (bmp->GetHeight())/10);
	  input_window->SetBmp(bmp);
	  output_window->SetBmp(NULL);
	  output_window->SetScrollbars(10,10,0,0);
	  input_window->Refresh();
	  output_window->Refresh();
	  // enable input fields
	  regions->Enable();
	  select_region_button->SetLabel("Select classes");  // reset button
							     // label
	  if (regs != NULL) 
	    {
	      delete [] regs;  // remove all rectangle selections
	      regs = NULL;
	      act_region = -1;
	      imageop->SetNoRegions(-1);
	    }
	  doit_button->Disable();
	  Refresh();
	}
    }
}


void MyFrame::OnSave(wxCommandEvent& event)
{
  if (imageop->IsOutput()) // if there is anything to save
    {
      wxString image_name;
      wxFileDialog* fdialog = new wxFileDialog(this, "Save file as", "", "", 
					       "BMP files (*.bmp)|*.bmp", 
					       wxSAVE | wxCHANGE_DIR | 
					       wxOVERWRITE_PROMPT);

      if (fdialog->ShowModal() == wxID_OK)
	{
	  image_name = fdialog->GetPath();
	  if (!imageop->SaveBmp(image_name)) // saving failed
	    wxMessageBox("Can't save image!", "ERROR");
	}
    }
}


void MyFrame::OnDoit(wxCommandEvent& event)
{
  wxString beta, t, T0, c, alpha;

  if ((beta=tbeta->GetValue()).Length() == 0)	
    {
      wxMessageBox("ß value missing", "Error");
      return;
    }
  else	// TODO: check value!
    imageop->SetBeta(atof(beta));
  if ((t=tt->GetValue()).Length() == 0)	
    {
      wxMessageBox("t value missing", "Error");
      return;
    }
  else	// TODO: check value!
    imageop->SetT(atof(t));
  if (op_choice->GetStringSelection() != "ICM")
    {
      if ((T0=tT0->GetValue()).Length() == 0)	
	{
	  wxMessageBox("T0 value missing", "Error");
	  return;
	}
      else	// TODO: check value!
	imageop->SetT0(atof(T0));
      if ((c=tc->GetValue()).Length() == 0)	
	{
	  wxMessageBox("c value missing", "Error");
	  return;
	}
      else	// TODO: check value!
	imageop->SetC(atof(c));
    }
  if (op_choice->GetStringSelection() == "MMD")
    {
      if ((alpha=talpha->GetValue()).Length() == 0)	
	{
	  wxMessageBox("alpha value missing", "Error");
	  return;
	}
      else	// TODO: check value!
	imageop->SetAlpha(atof(alpha));
    }

  timer_valid = FALSE; // timer's value is invalid. Used by GetTimer()
  Refresh();
  timer.Reset();       // reset timer
  timer.Start();       // start timer
  if (op_choice->GetStringSelection() == "Metropolis")
    {
      imageop->Metropolis();
    }
  else if (op_choice->GetStringSelection() == "MMD")
    {
      imageop->Metropolis(true);
    }
  else if (op_choice->GetStringSelection() == "ICM")
    {
      imageop->ICM();
    }
  else if (op_choice->GetStringSelection() == "Gibbs sampler")
    {
      imageop->Gibbs();
    }
  timer.Stop();       // stop timer
  timer_valid = TRUE; // timer's value is valid. Used by GetTimer()
  Refresh();
}


/* selection list handler
 */
void MyFrame::OnChoice(wxCommandEvent& event)
{
	if (op_choice->GetStringSelection() == "ICM")
	{
		tT0->Hide();
		tc->Hide();
	}
	else
	{
		tT0->Show();
		tc->Show();
	}
	if (op_choice->GetStringSelection() == "MMD")
	{
		talpha->Show();
	}
	else
	{
		talpha->Hide();
	}
	Refresh(); 
}


/* called whenever number of regions has been changed.
 */
void MyFrame::OnRegions(wxCommandEvent& event)
{
  select_region_button->SetLabel("Select classes"); // reset button label
  {
    delete [] regs;	// remove all rectangle selections
    regs = NULL;
    act_region = -1;
    imageop->SetNoRegions(-1);
    gaussians->SetValue("#\tMean\t\tVariance\n"); // clear parameter textfield

  }
  doit_button->Disable();
  wxString str = regions->GetValue(); // get entered value
  // TODO: check wheter the value is a positive integer!
  if (str.Length() != 0) select_region_button->Enable();
  else select_region_button->Disable();

}

		
/* on  "Select classes" button pushed
 */
void MyFrame::OnSelectRegion(wxCommandEvent& event)
{
  static int no_regions=0;
  if (select_region_button->GetLabel() == "Select classes")
    {
      act_region = 0;
      no_regions = atoi(regions->GetValue());     // get number of regions
      imageop->SetNoRegions(no_regions);
      regs = new int[no_regions*4];                   // aloccate memory
      for (int i=0; i<no_regions*4; ++i) regs[i] = 0; // init with 0
      select_region_button->SetLabel(act_region == no_regions-2? 
				     "Finish":"Next class");  // change label
    }
  else	// select next training rectangle
    {
      if (act_region == no_regions-2)
	select_region_button->SetLabel("Finish");  // change label
      if (IsSelected(act_region)) 
	imageop->CalculateMeanAndVariance(act_region);
      if (++act_region >= imageop->GetNoRegions())   // no more regions
	act_region = -1;                             // no act_region
      if (AllRegionsSelected()) {                    
	doit_button->Enable();            // Enable DoIt button
	select_region_button->Disable();  // Disable "Next class" button
      }
      
    }
  input_window->Refresh();
}


/* return the coordinates of the given training rectangle 
 */
void MyFrame::GetRegion(int &x, int &y, int &w, int &h, int region)
{
  x = y = w = h = 0;
  if (region == -1) // -1 ==> return current rectangle
    region = act_region;	
  if (region != -1) // otherwise compute coordinates
    {
      int x1 = regs[region*4];
      int y1 = regs[region*4+1];
      int x2 = regs[region*4+2];
      int y2 = regs[region*4+3];
      x = x1<x2 ? x1 : x2;		// left X coordinate
      y = y1<y2 ? y1 : y2;		// lower Y coordinate
      w = abs(x1-x2);			// width
      h = abs(y1-y2);			// height
    }
}


/*********************************************************************
/* Functions of ImageOperations class
/********************************************************************/
ImageOperations::ImageOperations(wxWindow *_frame)
{
  frame = _frame;
  in_image = out_image = NULL;

  no_regions = -1; // -1 ==> num. of regions has not been specified yet!
  beta = -1;
  T0 = -1;
  c = -1;
  K = 0;
  E = 0;
  T = 0;
  mean = variance = NULL;
  alpha = 0.1;
}


wxImage *ImageOperations::LoadBmp(wxString bmp_name)
{
  wxImage *img = new wxImage(bmp_name);
  if (img->Ok()) // set new values						
    {
      in_image = img;
      height = in_image->GetHeight();
      width = in_image->GetWidth();
      out_image = NULL;
    }
  return in_image;
}


bool ImageOperations::SaveBmp(wxString bmp_name)
{
  return out_image->SaveFile(bmp_name, wxBITMAP_TYPE_BMP);
}


bool ImageOperations::IsOutput()
{
  if (out_image == NULL) return false;
  else return true;
}


void ImageOperations::SetNoRegions(int n) 
{ 
  no_regions = n;
  if (n == -1)
    {
      delete [] mean;
      delete [] variance;
      mean = variance = NULL;
    }
  else 
    {
      mean = new double[n]; 
      variance = new double[n]; 
      for (int i=0; i<n; ++i) mean[i] = variance[i] = -1;
    }
}

	
/* Compute mean and variance for a given region
 */
void ImageOperations::CalculateMeanAndVariance(int region)
{
  if (in_image != NULL)
    {
      int x, y, w, h;
      int i, j;
      ((MyFrame *)frame)->GetRegion(x, y, w, h, region);

      double sum = 0, sum2=0;
      unsigned char *in_data = in_image->GetData();
      for (i=y; i<y+h; ++i)
	for (j=x; j<x+w; ++j) {
	  sum += in_data[i*width*3+j*3];
	  sum2 += in_data[i*width*3+j*3]*in_data[i*width*3+j*3];
	}
      mean[region] = sum/(w*h);
      variance[region] = (sum2 - (sum*sum)/(w*h))/(w*h-1);
//       sum = 0;
//       for (i=y; i<y+h; ++i)
// 	for (j=x; j<x+w; ++j)
// 	  sum += pow(in_data[i*width*3+j*3] - mean[region], 2);
//       variance[region] = sum/(w*h);
      if (variance[region] == 0) variance[region] = 1e-10;
      // print parameters in gaussians textfield
      *gaussians << region+1 << "\t" << mean[region] << "\t\t" << variance[region] << "\n";
    }
}


double ImageOperations::Singleton(int i, int j, int label)
{
  return log(sqrt(2.0*3.141592653589793*variance[label])) +
    pow((double)in_image_data[i][j]-mean[label],2)/(2.0*variance[label]);
}


double ImageOperations::Doubleton(int i, int j, int label)
{
  double energy = 0.0;

  if (i!=height-1) // south
    {
      if (label == classes[i+1][j]) energy -= beta;
      else energy += beta;
    }
  if (j!=width-1) // east
    {
      if (label == classes[i][j+1]) energy -= beta;
      else energy += beta;
    }
  if (i!=0) // nord
    {
      if (label == classes[i-1][j]) energy -= beta;
      else energy += beta;
    }
  if (j!=0) // west
    {
      if (label == classes[i][j-1]) energy -= beta;
      else energy += beta;
    }
  return energy;
}


/* compute global energy
 */
double ImageOperations::CalculateEnergy()
{
  double singletons = 0.0;
  double doubletons = 0.0;
  int i, j, k;
  for (i=0; i<height; ++i)
    for (j=0; j<width; ++j)
      {
	k = classes[i][j];
	// singleton
	singletons += Singleton(i,j,k);
	// doubleton
	doubletons += Doubleton(i,j,k); // Note: here each doubleton is
					// counted twice ==> divide by
					// 2 at the end!
      }
  return singletons + doubletons/2; 
}



double ImageOperations::LocalEnergy(int i, int j, int label)
{
  return Singleton(i,j,label) + Doubleton(i,j,label);
}


/* Initialize segmentation
 */
void ImageOperations::InitOutImage()
{
  int i, j, r;
  double e, e2;	 // store local energy

  unsigned char *in_data = in_image->GetData();

  in_image_data = new int* [height]; // allocate memory for in_image_data
  for (i=0; i<height; ++i)
    in_image_data[i] = new int[width];

  for (i=0; i<height; ++i)
    for (j=0; j<width; ++j)
      in_image_data[i][j] = in_data[(i*width*3)+j*3];

	
  classes = new int* [height]; // allocate memory for classes
  for (i=0; i<height; ++i)
    classes[i] = new int[width];
  /* initialize using Maximum Likelihood (~ max. of singleton energy)
   */
  for (i=0; i<height; ++i)
    for (j=0; j<width; ++j)
      {
	e = Singleton(i, j, 0);
	classes[i][j] = 0;
	for (r=1; r<no_regions; ++r)
	  if ((e2=Singleton(i, j, r)) < e)
	    {
	      e = e2;
	      classes[i][j] = r;
	    }
      }
}

/* Create and display the output image based on the current labeling.
 * Executed at each iteration.
 */
void ImageOperations::CreateOutput()
{
  int i, j;
  unsigned char *out_data;

  /* Do not count GUI overhead
   */
  timer.Stop();

  out_data = (unsigned char *)malloc(width*height*3*sizeof(unsigned char));
	
  // cretae output image
  for (i=0; i<height; ++i)
    for (j=0; j<width; ++j)
      {
	out_data[(i*width*3) + j*3] = 
	  (unsigned char)(classes[i][j]*255/no_regions);
	out_data[(i*width*3) + j*3 +1] = 
	  (unsigned char)(classes[i][j]*255/no_regions);
	out_data[(i*width*3) + j*3 +2] = 
	  (unsigned char)(classes[i][j]*255/no_regions);
      }

  free (out_image);
  out_image = new wxImage(width, height, out_data);

  // and display it
  ((MyFrame *)frame)->GetOutputWindow()->SetScrollbars(10,10,(out_image->GetWidth())/10,(out_image->GetHeight())/10);
  ((MyFrame *)frame)->GetOutputWindow()->SetBmp(out_image);
  ((MyFrame *)frame)->GetOutputWindow()->Refresh();
  ((MyFrame *)frame)->GetOutputWindow()->Update();
  frame->RefreshRect(wxRect(645, 360, 100, 100));
  frame->Update();
  /* Continue timer
   */
  timer.Start();
}


/* Metropolis & MMD
 */
void ImageOperations::Metropolis(bool mmd)
{
  InitOutImage();
  int i, j;
  int r;
  double kszi = log(alpha);  // This is for MMD. When executing
			    // Metropolis, kszi will be randomly generated.
  double summa_deltaE;
  
  TRandomMersenne rg(time(0));  // create instance of random number generator

  K = 0;
  T = T0;
  E_old = CalculateEnergy();

  do
    {
      summa_deltaE = 0.0;
      for (i=0; i<height; ++i)
	for (j=0; j<width; ++j)
	  {
	    /* Generate a new label different from the current one with
	     * uniform distribution.
	     */
	    if (no_regions == 2)
	      r = 1 - classes[i][j];
	    else
	      r = (classes[i][j] +
		   (int)(rg.Random()*(no_regions-1))+1) % no_regions;
	    if (!mmd)  // Metropolis: kszi is a  uniform random number
	      kszi = log(rg.Random()); 
	    /* Accept the new label according to Metropolis dynamics.
	     */
	    if (kszi <= (LocalEnergy(i, j, classes[i][j]) -
			 LocalEnergy(i, j, r)) / T) {
	      summa_deltaE += 
		fabs(LocalEnergy(i, j, r) - LocalEnergy(i, j, classes[i][j]));
	      E_old = E = E_old - 
		LocalEnergy(i, j, classes[i][j]) + LocalEnergy(i, j, r);
	      classes[i][j] = r;
	    }
	  }
      T *= c;         // decrease temperature
      ++K;	      // advance iteration counter
      CreateOutput(); // display current labeling
    } while (summa_deltaE > t); // stop when energy change is small
}


/* ICM
 */
void ImageOperations::ICM()
{
  InitOutImage();
  int i, j;
  int r;
  double summa_deltaE;

  K = 0;
  E_old = CalculateEnergy();

  do
    {
      summa_deltaE = 0.0;
      for (i=0; i<height; ++i)
	for (j=0; j<width; ++j)
	  {
	    for (r=0; r<no_regions; ++r)
	      {
		if (LocalEnergy(i, j, classes[i][j]) > LocalEnergy(i, j, r))
		  {
		    classes[i][j] = r;
		  }
	      }
	  }
      E = CalculateEnergy();
      summa_deltaE += fabs(E_old-E);
      E_old = E;

      ++K;	      // advance iteration counter
      CreateOutput(); // display current labeling
    }while (summa_deltaE > t); // stop when energy change is small
}


/* Gibbs sampler
 */
void ImageOperations::Gibbs()
{
  InitOutImage();
  int i, j;
  double *Ek;		       // array to store local energies
  int s;
  double summa_deltaE;
  double sumE;
  double z;
  double r;

  TRandomMersenne rg(time(0)); // make instance of random number generator

  Ek = new double[no_regions];

  K = 0;
  T = T0;
  E_old = CalculateEnergy();

  do
    {
      summa_deltaE = 0.0;
      for (i=0; i<height; ++i)
	for (j=0; j<width; ++j)
	  {
	    sumE = 0.0;
	    for (s=0; s<no_regions; ++s)
	      {
		Ek[s] = exp(-LocalEnergy(i, j, s)/T);
		sumE += Ek[s];
	      }
	    r = rg.Random();	// r is a uniform random number
	    z = 0.0;
	    for (s=0; s<no_regions; ++s)
	      {
		z += Ek[s]/sumE; 
		if (z > r) // choose new label with probabilty exp(-U/T).
		  {
		    classes[i][j] = s;
		    break;
		  }
	      }
	  }
      E = CalculateEnergy();
      summa_deltaE += fabs(E_old-E);
      E_old = E;

      T *= c;         // decrease temperature
      ++K;	      // advance iteration counter
      CreateOutput(); // display current labeling
    } while (summa_deltaE > t); // stop when energy change is small

  delete Ek;
}
