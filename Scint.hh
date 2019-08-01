#ifndef __SCINT__
#define __SCINT__ 1

#include <iostream>
#include <vector>
#include <string>

#include <TVector3.h>
#include <TMath.h>
#include <TH1D.h>

class Line
{
public:
  std::string name;
  TVector3 slope, point;
  double kParall=-999999999.;
  
  Line(TVector3 u, TVector3 p, std::string n=""):name(n),slope(u),point(p)
  {}
  
  TVector3 operator()(double t)
  {
    return slope * t + point;
  }

  TVector3 Intersection(Line l)
  {
    double den = slope.Cross(l.slope).Mag2();
    double s = den > 0.?((l.point - point).Cross(l.slope)) * (slope.Cross(l.slope)) / den : kParall;
    return point + slope * s;
  }
};

class Surface
{
public:
  std::string name;
  TVector3 normal, point;
  
  Surface(TVector3 n, TVector3 p, std::string m=""):name(m),normal(n),point(p)
  {}
  Surface(TVector3 n, TVector3 p, const char* m):name(m),normal(n),point(p)
  {}

  TVector3 operator()(double, double)
  {
    std::cerr<<"Surface operator() NOT IMPLEMENTED"<<std::endl;
    TVector3 null(0.0,0.0,0.0);
    return null;
  }

  TVector3 Intersection(Line l)
  {
    double den = l.slope.Dot(normal);
    double t = (den != 0.) ? ((point - l.point) * normal / den) : l.kParall;
    return l(t);
  }

  Line Intersection(Surface p)
  {
    std::cerr<<"Surface Intersection(Surface ) NOT IMPLEMENTED"<<std::endl;
    TVector3 null(0.0,0.0,0.0);
    Line zero(null,null);
    return zero;
  }
};

class SiPM
{
public:
  double posz;
  double sizex;
  double sizey;
  double halfsizex;
  double halfsizey;
  std::vector<Surface> surfaces;
  TVector3 hit;

  SiPM(double z, double x, double y): posz(z), sizex(x), sizey(y)
  {
    surfaces.emplace_back(TVector3(0.,0.,-1.),TVector3(0.,0.,posz),"SiPMa");
    surfaces.emplace_back(TVector3(0.,0.,1.),TVector3(0.,0.,-posz),"SiPMb");
    halfsizex=0.5*sizex;
    halfsizey=0.5*sizey;
  }

  bool Hit(double x, double y, double z)
  {
    if( (fabs(x) < halfsizex) && (fabs(y) < halfsizey) && (fabs(z) >= posz) )
      {
	hit.SetXYZ(x,y,z);
	return true;
      }
    else return false;
  }

  bool Hit(TVector3 p)
  {
    return Hit(p.X(),p.Y(),p.Z());
  }

  bool Hit(Line photon)
  {
    return Hit(surfaces[0].Intersection( photon )) || Hit(surfaces[1].Intersection( photon ));
  }
};

class Scint
{
public:
  double length; // z [mm]
  double sidex;
  double sidey;

private:
  std::vector<Surface> surfaces;
  //double lim[4];

  SiPM det;

  TH1D* h_internal_reflection;

  Line null;

public:
  double index_of_refraction;
  double speed_of_light;

  double reflectivity;
  
  // ALPHA-g geometry of the bar
  Scint():length(2500.),sidex(20.),sidey(22.),// mm
	  det(length*0.5,sidex,sidey),
	  null(TVector3(0.,0.,0.),TVector3(0.,0.,0.),"null"),
	  index_of_refraction(1.58), reflectivity(0.99)
  {
    speed_of_light = TMath::C()*1.e-06 / index_of_refraction; // mm/ns
    Surfaces();
    h_internal_reflection = new TH1D("h_internal_reflection",
				     "Internal Reflection Angle;#theta_{i} [deg]",
				     1000,0.,90.);
  }
  
  void Surfaces()
  {
    double half_length = 0.5*length,
      halfsidex = 0.5*sidex, halfsidey = 0.5*sidey;
    
    surfaces.emplace_back(TVector3(0.,1.,0.),
			  TVector3(0.,-halfsidey,half_length),
			  "x-z plane @ -y");

    surfaces.emplace_back(TVector3(-1.,0.,0.),
			  TVector3(halfsidex,0.,half_length),
			  "y-z plane @ +x");

    surfaces.emplace_back(TVector3(0.,-1.,0.),
			  TVector3(0.,halfsidey,half_length),
			  "x-z plane @ +y");

    surfaces.emplace_back(TVector3(1.,0.,0.),
			  TVector3(-halfsidex,0.,half_length),
			  "y-z plane @ -x");

    //    double lim[4] = {halfsidex,halfsidey,halfsidex,halfsidey};
    //for(int i=0;i<4;++i) lim[i]=((i%2)==0)?halfsidex:halfsidey;
  }

  std::vector<TVector3> Propagate(Line init);
  Line Reflection(Line photon);
  bool Inside(bool yplane, TVector3 r)
  {
    if( yplane && (TMath::Abs( r.X() ) < 0.5*sidex) ) return true;
    else if( !yplane && (TMath::Abs( r.Y() ) < 0.5*sidey) ) return true;
    else return false;
  }
  
  void SmartPrint(TVector3 v, TString info)
  {
    std::cout<<info<<" = ("<<v.X()<<","<<v.Y()<<","<<v.Z()<<")"<<std::endl;
  }

  TH1D* GetReflAngleHisto() { return h_internal_reflection; }
};
#endif
