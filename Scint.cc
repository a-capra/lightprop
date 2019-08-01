#include "Scint.hh"
extern bool bounce_info;
#include <TRandom3.h>
extern TRandom3 rndm;

#include <cassert>

Line Scint::Reflection(Line photon)
{
  if( rndm.Uniform() > reflectivity ) return null;
  TVector3 intersect;
  int s=-1;
  bool yplane;
  for(int i=0;i<4;++i)
    { 
      // bounce
      intersect = surfaces[i].Intersection( photon );
      
      // avoids numerical errors
      if( (intersect - photon.point).Mag() < 1.e-6 ) continue;

      // photon and surface are parallel
      if( intersect == photon(photon.kParall) ) continue;

      // ensure propagation (forward)
      if( photon.slope.Z() > 0. && intersect.Z() < photon.point.Z() ) continue;
      if( photon.slope.Z() < 0. && intersect.Z() > photon.point.Z() ) continue;

      // indetify bouncing surface
      yplane = ( (i%2) == 0 );
      if( Inside( yplane, intersect) ) 
	{
	  s=i;
      	  break;
	}
    }
  
  if( s < 0 ) 
    {
      std::cerr<<" Reflection fail"<<std::endl;
      return null;
    }
  //assert(!(s<0));
  
  TVector3 uprime;
  if( yplane )
    uprime.SetXYZ(photon.slope.X(),-1.*photon.slope.Y(),photon.slope.Z());
  else
    uprime.SetXYZ(-1.*photon.slope.X(),photon.slope.Y(),photon.slope.Z());

  //double thetai = (surfaces[s].normal.Angle(uprime) - TMath::PiOver2())*TMath::RadToDeg();
  double thetai = surfaces[s].normal.Angle(uprime) * TMath::RadToDeg();
  h_internal_reflection->Fill( thetai );
  if( bounce_info ) std::cout<<"Incidence angle "<<thetai<<" deg"<<std::endl;
  Line reflray(uprime,intersect,photon.name);
  return reflray;
}

std::vector<TVector3> Scint::Propagate(Line ray)
{
  std::vector<TVector3> refelections;
  int N=0;
  while(1)
    {
      if( bounce_info )
	{
	  std::cout<<"***** Bounce # "<<N<<" for "<<ray.name<<" *****"<<std::endl;
	  SmartPrint(ray.slope,"direction");
	  SmartPrint(ray.point,"position");
	}
      // if( N > 100000 )
      // 	{
      // 	  std::cout<<"KILL ME"<<std::endl;
      // 	  refelections.clear();
      // 	  break;
      // 	}
      refelections.push_back( ray.point );
      if( det.Hit( ray ) )
	{
	  refelections.push_back( det.hit );
	  if( bounce_info ) SmartPrint(det.hit,"HIT");
	  break;
	}
      else
	ray = Reflection( ray );
      if( bounce_info ) std::cout<<"******************"<<std::endl;
      if( ray.slope.Z() == null.slope.Z() ) 
	{
	  if( bounce_info ) std::cout<<"Absorbed after "<<N<<" reflections\n"<<std::endl;
	  refelections.clear();
	  break;
	}
      ++N;
    }
  return refelections;
}
