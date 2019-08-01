#include "Scint.hh"

#include <TRandom3.h>
TRandom3 rndm(19850618);
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TPolyLine3D.h>
//#include <TBRIK.h>
#include <TAxis3D.h>

using namespace TMath;

bool bounce_info=false;

using namespace std;
int main(int argc, char** argv)
{
  int Nph = 100000; // number of photons
double z0=0.;
  if( argc == 2 )
    {
      Nph = atoi(argv[1]);
    }
  else if( argc == 3 )
    {
      Nph = atoi(argv[1]);
z0 = (double) atof(argv[2]);
    }
  else if( argc == 4 )
    {
      Nph = atoi(argv[1]);
z0 = (double) atof(argv[2]);
      bounce_info = (bool) atoi(argv[3]);
    }

  // my scintillating bar
  Scint bar;
  
  // photon
  double ux,uy,uz; // slope
double x0=0., y0=0.;//,z0;// point

  TString fname  = TString::Format("LTT_barZ%1.f_X%1.f_Y%1.f_z0%1.1f_n%1.2f_refl%0.2f.root",
				     bar.length,bar.sidex,bar.sidey,z0,
				     bar.index_of_refraction,
				     bar.reflectivity);

  TFile* fout = TFile::Open(fname.Data(),"RECREATE");
  
  TH1D* ht = new TH1D("ht","Light Transit Time;t [ns];#gamma",20000,0.,200.);
  TH2D* hxy = new TH2D("hxy","Light Position on SiPM;x [mm];y [mm];#gamma",
		       100,-0.5*bar.sidex,0.5*bar.sidex,100,-0.5*bar.sidey,0.5*bar.sidey);
  TH1D* hNr = new TH1D("hNr","Number of Reflections",2000,0.,2000.);
  TH1D* htheta = new TH1D("htheta","Photon Angle w.r.t. Bar Axis;#theta [rad];#gamma",1000,
			  0.,Pi());
  TH2D* hthetarefl = new TH2D("hthetarefl",
			      "Photon Angle w.r.t. Bar Axis Vs # of Reflections;#theta [rad];Number of Reflections;#gamma",
			      100,0.,Pi(),
			      2000,0.,2000.);
  TH1D* hdist = new TH1D("hdist","Light Travelled Distance;d [mm];#gamma",1000,0.,1.e4);
  
  TH2D* hzt = new TH2D("hzt","Time on each end; TOP/BOT;t [ns]",2,-1.,1.,20000,0.,200.);
  //TH2D* hzt = new TH2D("hzt","Time on each end; TOP/BOT;t [ns]",3,-1300.,1300.,20000,0.,200.);

  TH1D* ht_bot = new TH1D("ht_bot","Light Transit Time Bottom;t [ns];#gamma",20000,0.,200.);
  TH2D* hxy_bot = new TH2D("hxy_bot","Light Position on Bottom SiPM;x [mm];y [mm];#gamma",
		       100,-0.5*bar.sidex,0.5*bar.sidex,100,-0.5*bar.sidey,0.5*bar.sidey);
  TH1D* ht_top = new TH1D("ht_top","Light Transit Time Top;t [ns];#gamma",20000,0.,200.);
  TH2D* hxy_top = new TH2D("hxy_top","Light Position on Top SiPM;x [mm];y [mm];#gamma",
		       100,-0.5*bar.sidex,0.5*bar.sidex,100,-0.5*bar.sidey,0.5*bar.sidey);

  std::vector<TPolyLine3D*> rays;
  // number of propagated photons
  int iph = 0;
  while(iph<Nph)
    {
      if( bounce_info || iph%100==0 )
	cout<<"photon # "<<iph<<endl;

      // generate direction
      rndm.Sphere(ux,uy,uz,1.);
      TVector3 u(ux,uy,uz);
      // generate position on bar
      //z0 = rndm.Uniform(-0.5*bar.length,0.5*bar.length);
      //z0=0.;
      TVector3 r0(x0,y0,z0);

if( bounce_info )
  cout<<"direction theta: "<<u.Theta()*RadToDeg()<<" deg"<<endl;
     
      string lname="photon ";
      lname+=to_string(iph);
      Line photon(u,r0,lname);
      vector<TVector3> refl = bar.Propagate( photon );
      
      if( refl.size() == 0 ) continue;

      TPolyLine3D* ray = new TPolyLine3D;
      int np=0;
      for(auto it=refl.begin(); it != refl.end(); ++it ) ray->SetPoint(np++,it->X(),it->Y(),it->Z());
rays.push_back(ray);

      double Nr=double(refl.size()-2),// number of reflections 
	dist = 0.; // mm
      hNr->Fill( Nr );
      for(size_t i=1; i<refl.size(); ++i)
	dist += (refl[i-1]-refl[i]).Mag();
      hdist->Fill(dist);
      double time = dist / bar.speed_of_light;
      ht->Fill( time );
      
      TVector3 f = refl.back();
      hxy->Fill( f.X(), f.Y() );

      htheta->Fill( u.Theta() );
      hthetarefl->Fill( u.Theta(), Nr );

      hzt->Fill( f.Z()/bar.length, time );
      //hzt->Fill( f.Z(), time );

      if( (f.Z() + 0.5*bar.length) < 1.e-9 )
	{
	  ht_bot->Fill( time );
	  hxy_bot->Fill( f.X(), f.Y() );
	}
      else if( (f.Z() - 0.5*bar.length) < 1.e-9 )
	{
	  ht_top->Fill( time );
	  hxy_top->Fill( f.X(), f.Y() );
	}
      else
	cerr<<"Error ph: "<<iph<<" hit: "<<f.Z()<<" mm"<<endl;
      
      if( bounce_info || iph%100==0 )      
	cout<<"========================================================================\n"<<endl;
     ++iph;
    }

  bar.GetReflAngleHisto()->Write();
  fout->Write();
  fout->Close();

if( rays.size() == 0 ) return 1;

 //  TApplication app("lightprop",&argc,argv);
//   TAxis3D rulers;
//   TCanvas* c1 = new TCanvas("c1","c1",1200,1000);
//   c1->cd();
//   //  rulers.Draw("ogl");
// rulers.Draw();
//   //  bar->Draw("oglsame");
//   //ray->Draw("oglsame");
//   for(auto it = rays.begin(); it != rays.end(); ++it)
//     (*it)->Draw("same");
//   app.Run();

  return 0;
}
