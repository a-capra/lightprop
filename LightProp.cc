//======================================================================================
// g++ -Wall -O3 `root-config --cflags` `root-config --glibs` -o LightProp LightProp.cc
//
//
//======================================================================================
#include <TMath.h>
using namespace TMath;
#include <TVector3.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TFile.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TPolyLine3D.h>
#include <TBRIK.h>
#include <TAxis3D.h>

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

// generate direction
//TRandom3 rndm(20160903);
TRandom3 rndm(19850618);

void SmartPrint(TVector3 v, TString info)
{
  cout<<info<<" = ("<<v.X()<<","<<v.Y()<<","<<v.Z()<<")"<<endl;
}

TVector3 Line(double t, TVector3 u, TVector3 r0)
{
  return u * t + r0;
}

TVector3 LinePlaneIntersection(TVector3 u, TVector3 r0,
			       TVector3 n, TVector3 p0)
{
  double t = -999999999.;
  if( u.Dot(n) != 0. )
    t = ( ( p0 - r0 ).Dot(n) ) / u.Dot(n);

  return Line(t,u,r0);
}

TH1D* h_internal_reflection;
bool bounce_info=false;
double reflectivity = 0.99;
// SiPM plane
TVector3 n0(0.,0.,-1.), p0;
//std::vector<TPolyLine3D*> rays;
double reflection(TVector3& SiPM, // output
		  double& number_of_reflections,// output
		  TVector3 u, TVector3 r0, // line vector
		  double s1, double s2, // half lengths
		  double length)
{
  double half_length = 0.5*length;

  TVector3 n[4], p[4];
  n[0].SetXYZ(0.,1.,0.); p[0].SetXYZ(0.,-s2,half_length); // x-z plane @ -y
  n[1].SetXYZ(-1.,0.,0.); p[1].SetXYZ(s1,0.,half_length); // y-z plane @ +x
  n[2].SetXYZ(0.,-1.,0.); p[2].SetXYZ(0.,s2,half_length); // x-z plane @ +y
  n[3].SetXYZ(1.,0.,0.); p[3].SetXYZ(-s1,0.,half_length); // y-z plane @ -x
  double lim[4] = {s1,s2,s1,s2};

  double distance=0., z=0.;
  int Nrefl=1;
  TVector3 r0prime=r0, uprime = u;
  // //  ray = new TPolyLine3D;
  // TPolyLine3D* ray = new TPolyLine3D;
  // ray->SetPoint(0,r0prime.X(),r0prime.Y(),r0prime.Z());
  // //rays.emplace_back();
  // //rays.back()->SetPoint(0,r0prime.X(),r0prime.Y(),r0prime.Z());
  // //rays.back()->SetLineColor(kBlue);
  // //rays.back()->SetLineWidth(2);

  while( z < length )
    {
      if( bounce_info )
	{
	  cout<<"***** Bounce *****"<<endl;
	  SmartPrint(uprime,"direction u\'");
	  SmartPrint(r0prime,"point r0\'");
	}

      int m=-1,mm=0;
      // calculate intersection of the light ray with
      // bar surface, a.k.a. photon bouncing point 
      for(int i=0; i<4; ++i)
	{
	  TVector3 r = LinePlaneIntersection(uprime, r0prime, n[i], p[i]);
	  //SmartPrint(r,"possible plane intersect");

	  if( (r - r0prime).Mag() < 1.e-6 ) // avoids numerical errors
	    continue;

	  TVector3 parall = Line(-999999999.,uprime,r0prime);
	  if( r == parall )
	    continue;

	  double comp;
	  if( (i%2) == 0 )
	    comp = Abs( r.X() );
	  else 
	    comp = Abs( r.Y() );

	  // of course only one point exists
	  // and is inside the bar, which is all at z>z_starting_point
	  if( comp < lim[i] && r.Z() > r0.Z() )
	    {
	       m = i;
	       ++mm;
	    }
	}

      if( mm > 1 ) // fail
	{
	  cerr<<"Error: mutliple intersections ("<<mm<<")"<<endl;
	  return -6.;
	}
      
      if( m < 0 ) // fail
	{
	  cerr<<"Error: no intersection ("<<m<<")"<<endl;
	  return -1.;
	}
      
      if( rndm.Uniform() > reflectivity )
	{
	  cout<<" $$$$$$$ Absorption! $$$$$$$"<<endl;
	  return -5.;
	}

      // Reflection from surface
      TVector3 intersect = LinePlaneIntersection(uprime, r0prime, n[m], p[m]);
      if( bounce_info )
	{
	  TString def_plane_int = TString::Format("%d Reflection from surface %d ", Nrefl, m);
	  SmartPrint(intersect,def_plane_int);
	}

      // path length
      distance += ( intersect - r0prime ).Mag();

      // reflection point
      r0prime = intersect;
      // ray->SetPoint(Nrefl,r0prime.X(),r0prime.Y(),r0prime.Z());
      // //rays.back()->SetPoint(Nrefl,r0prime.X(),r0prime.Y(),r0prime.Z());

      // incident angle
      //double thetai = uprime.Angle(n[m]) - PiOver2();
      double thetai = n[m].Angle(uprime);
      if( bounce_info )
	cout<<"Incidence angle "<<thetai*RadToDeg()<<" deg ("<<(uprime.Angle(n[m])-PiOver2())*RadToDeg()<<" deg)"<<endl;
      
      // calculate direction of reflected ray
      // this is the right formula but it doesn't work
      // uprime = uprime + 2. * ( n[m].Dot(uprime) ) * n[m];
      // instead use this
      TVector3 utemp(uprime);
      if( (m%2) == 0 ) // normal along y
      	uprime.SetXYZ( utemp.X(), -1.*utemp.Y(), utemp.Z() );
      else
      	uprime.SetXYZ( -1.*utemp.X(), utemp.Y(), utemp.Z() ); 
      
      if( uprime.Z() <= 0. )// fail
	{
	  cerr<<"Error: Backward reflection v1"<<endl;
	  return -2.;
	}
      if( uprime.Dot(n0) >= 0. )// the same fail
	{
	  cerr<<"Error: Backward reflection v2"<<endl;
	  return -3.;
	}

      // reflected angle
      //double thetaf = n[m].Angle(-uprime) - PiOver2();
      double thetaf = n[m].Angle(-uprime);
      if( thetai != thetaf ) // worst fail
      	{
	  cerr<<"Error: reflected angle "<<thetaf*RadToDeg()<<" != incident angle "<<thetai*RadToDeg()<<" deg"<<endl;
      	  return -4.;
      	}

      h_internal_reflection->Fill( thetai*RadToDeg() );

      // intersection with SiPM plane
      TVector3 f = LinePlaneIntersection(uprime, r0prime, n0, p0);
      if( bounce_info )
	{
	  SmartPrint(uprime,"reflected ray direction u\'");
	  SmartPrint(r0prime,"            from point r0\'");
	  SmartPrint(f,TString("SiPM plane f"));
	}

      if( Abs( f.X() ) < s1 && Abs( f.Y() ) < s2 )
	{
	  // cout<<"\t\tit\'s IN!"<<endl;
	  SiPM = f;
	  distance += (f-r0prime).Mag();
	  z = f.Z();
	  //	  break;
	}
      else
	{
	  // cout<<"\t\tit\'s OUT..."<<endl;
	  z = intersect.Z();
	  ++Nrefl;
	}

      if( bounce_info ) cout<<"******************"<<endl;
    }
  //  rays.push_back(ray);
  
  number_of_reflections = double(Nrefl);
  return distance;
}

int main(int argc, char** argv)
{
  bool alpha2=false;
  int Nph = 100000; // number of photons
  if( argc == 2 )
    alpha2 = (bool) atoi(argv[1]);
  else if( argc == 3 )
    {
      alpha2 = (bool) atoi(argv[1]);
      Nph = atoi(argv[2]);
    }
  else if( argc == 4 )
    {
      alpha2 = (bool) atoi(argv[1]);
      Nph = atoi(argv[2]);
      bounce_info = (bool) atoi(argv[3]);
    }

  // geometry of the bar (ALPHA-g)
  double length = 2500., // mm
    sidex = 20., sidey=22.;

  // (ALPHA-2)
  if( alpha2 )
    {
      length = 1650.;
      sidex *= 2.;
    }
  double halfsidex = 0.5*sidex, halfsidey = 0.5*sidey;

  // SiPM position
  p0.SetXYZ(0.,0.,length);
  //  p0.SetXYZ(0.,0.,0.5*length);

  // line
  double ux,uy,uz; // slope
  double x0=0., y0=0., // point
    z0 = length*0.5; 
  TVector3 r0(x0,y0,z0);

  // outside the bar already!
  if( Abs( x0 ) < halfsidex && Abs( y0 ) < halfsidey )
    if( z0 < 0. || z0 > length ) return -1;
  
  double index_of_refraction = 1.58,
			 c = C()*1.e-06; // mm/ns
  double speed_of_light = c / index_of_refraction;

  // number of photons to be propagated
  int iph = 0;

  TString fname  = TString::Format("LTT_barZ%1.f_X%1.f_Y%1.f_phx%1.2f_phy%1.2f_phz%1.2f_n%1.2f_refl%0.2f.root",
				   length,sidex,sidey,
				   x0, y0, z0,
				   index_of_refraction,
				   reflectivity);

  TFile* fout = TFile::Open(fname.Data(),"RECREATE");
  TH1D* ht = new TH1D("ht","Light Transit Time;t [ns];#gamma",20000,0.,200.);
  TH1D* ht_direct = new TH1D("ht_direct","Direct Light Transit Time;t [ns];#gamma",
			     20000,0.,200.);
  TH2D* hxy = new TH2D("hxy","Light Position on SiPM;x [mm];y [mm];#gamma",
		       100,-halfsidex,halfsidex,100,halfsidey,halfsidey);
  TH2D* hxy_direct = new TH2D("hxy_direct","Direct Light Position on SiPM;x [mm];y [mm];#gamma",
			      100,-halfsidex,halfsidex,100,halfsidey,halfsidey);
  TH1D* hNr = new TH1D("hNr","Number of Reflections",2000,0.,2000.);
  TH1D* htheta = new TH1D("htheta","Photon Angle w.r.t. Bar Axis;#theta [rad];#gamma",1000,
			  0.,PiOver2());
  TH1D* htheta_direct = new TH1D("htheta_direct",
				 "Direct Photon Angle w.r.t. Bar Axis;#theta [rad];#gamma",
				 1000,0.,PiOver2());
  TH2D* hthetarefl = new TH2D("hthetarefl",
			      "Photon Angle w.r.t. Bar Axis Vs # of Reflections;#theta [rad];Number of Reflections;#gamma",
			      100,0.,PiOver2(),
			      2000,0.,2000.);
  TH1D* hdist = new TH1D("hdist","Light Travelled Distance;d [mm];#gamma",1000,0.,1.e4);

  h_internal_reflection = new TH1D("h_internal_reflection","Internal Reflection Angle;#theta_{i} [deg]",1000,0.,90.);
  

  // loop over N photons
  while( iph < Nph )
    {
      //      cout<<"photon # "<<iph;
      cout<<"photon # "<<iph<<endl;

      // generate direction -> line slope
      rndm.Sphere(ux,uy,uz,1.);
      TVector3 u(ux,uy,uz);
      
      TVector3 f;
      if( u.Dot(n0) < 0. )
	{
	  // project/propagate light to end
	  // line - plane intersection
	  f = LinePlaneIntersection(u,r0,n0,p0);
	}
      else
	{
	  cout<<" backwards\n";
	  continue;
	}

      // cout<<"\tits direction: ("<<ux<<","<<uy<<","<<uz<<")";
      
      TVector3 parall = Line(-999999999.,u,r0);
      if( f == parall )
	{
	  cout<<" parallel\n";
	  continue;
	}

      // cout<<"\tintersection ("<<f.X()<<","<<f.Y()<<","<<f.Z()<<")\tdirect? ";

      htheta->Fill( u.Theta() );

      double Nr=0.,// number of reflections
	time = 0., d; // ns
      TVector3 interSiPM; // interesection with SiPM plane
	
      // if intersect -> end
      if( Abs( f.X() ) < halfsidex && Abs( f.Y() ) < halfsidey )
	{
	  // cout<<"yes"<<endl;
	  d = (f-r0).Mag();
	  time = d / speed_of_light;
	  
	  ht_direct->Fill(time);
	  hxy_direct->Fill( f.X(), f.Y() );
	  htheta_direct->Fill( u.Theta() );

	  interSiPM = f;
	}
      else // else bounce -> reflection until end
	{
	  // cout<<"no"<<endl;
	  	  
	  //SmartPrint(u,"direction u");
	  d = reflection( interSiPM,
				 Nr,
				 u, r0,
				 halfsidex, halfsidey, length );
	  //cout<<" Total Distance: "<<d<<" mm"<<endl;
	  if( d < 0. )
	    {
	      //	      cout<<"\n";
	      continue;
	    }
	  cout<<" Total Distance: "<<d<<" mm"<<endl;
	  time = d / speed_of_light;
	}

      ht->Fill( time );
      hxy->Fill( interSiPM.X(), interSiPM.Y() );

      hNr->Fill( Nr );

      hthetarefl->Fill( u.Theta(), Nr );

      hdist->Fill(d);

      ++iph;
      //cout<<"\n"<<endl;
      cout<<"========================================================================\n"<<endl;
    }

  fout->Write();
  fout->Close();


  // //TBRIK* bar = new TBRIK("bar","bar",0,halfsidex,halfsidey,0.5*length);
  // // TBRIK* bar = new TBRIK("bar","bar",0,halfsidex,halfsidey,length);
  // // bar->SetFillStyle(0);
  // TApplication app("lightprop",&argc,argv);
  // TAxis3D rulers;
  // TCanvas* c1 = new TCanvas("c1","c1",1200,1000);
  // c1->cd();
  // //  rulers.Draw("ogl");
  // rulers.Draw();
  // //  bar->Draw("oglsame");
  // //ray->Draw("oglsame");
  // for(auto it = rays.begin(); it != rays.end(); ++it)
  //   (*it)->Draw("same");
  // app.Run();

  return 0;
}
