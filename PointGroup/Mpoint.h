#ifndef _Mpoint_h_
#define _Mpoint_h_

#include <stdlib.h>
#include <assert.h>
#include <Misc.h>

class Mpoint{
 public:
  Mpoint(real y=0, real x=0, real vy=0, real vx=0, int nump=1, int numc=1) {
    Uy=vy;
    Ux=vx;
    X=x;
    Y=y;
	NumPoints=nump;
	Vy=new real[numc];
	Vx=new real[numc];
	Prob=new real[numc];
	int i;
	for(i=0; i<NumPoints; ++i) {
		Vy[i]=0;
		Vx[i]=0;
		Prob[i]=1.0/numc;
	}
    Link=new real[NumPoints];
	for(i=0; i<NumPoints; ++i)
		Link[i]=1.0;
  }

  Mpoint(const Mpoint& cp) {
    Ux=cp.GetVelocityY();
    Uy=cp.GetVelocityX();
    X=cp.GetX();
    Y=cp.GetY();

	NumCandidates=cp.GetNumCandidates();
	Vy=new real[NumCandidates];
	Vx=new real[NumCandidates];
	Prob=new real[NumCandidates];
	int i;
	for(i=0; i<NumCandidates; ++i) {
		Vy[i]=cp.GetEstimateVelocityY(i);
		Vx[i]=cp.GetEstimateVelocityX(i);
		Prob[i]=cp.GetProb(i);
	}

    NumPoints=cp.GetNumPoints();
	if(Link)
		delete [] Link;
	Link=new real[NumPoints];
 	for(i=0; i<NumPoints; ++i)
	   Link[i]=cp.GetLink(i);
  }

  const Mpoint& operator=(const Mpoint& cp) {
    Ux=cp.GetVelocityY();
    Uy=cp.GetVelocityX();
    X=cp.GetX();
    Y=cp.GetY();

	NumCandidates=cp.GetNumCandidates();
	Vy=new real[NumCandidates];
	Vx=new real[NumCandidates];
	Prob=new real[NumCandidates];
	int i;
	for(i=0; i<NumCandidates; ++i) {
		Vy[i]=cp.GetEstimateVelocityY(i);
		Vx[i]=cp.GetEstimateVelocityX(i);
		Prob[i]=cp.GetProb(i);
	}

    NumPoints=cp.GetNumPoints();
	if(Link)
		delete [] Link;
	Link=new real[NumPoints];
 	for(i=0; i<NumPoints; ++i)
	   Link[i]=cp.GetLink(i);

    return *this;
  }

  inline real GetVelocityY() const {return Uy;}
  inline real GetVelocityX() const {return Ux;}
  inline real GetX() const {return X;}
  inline real GetY() const {return Y;}
  inline int GetNumPoints() const {return NumPoints;}
  inline int GetNumCandidates() const {return NumCandidates;}
  inline real GetLink(int i) const {return Link[i];}
  inline real GetProb(int i) const {return Prob[i];}
  inline real GetEstimateVelocityY(int i) const {return Vy[i];}
  inline real GetEstimateVelocityX(int i) const {return Vx[i];}

  inline void SetEstimateVelocityY(int n, real v) {
    Vy[n]=v;
  }
  inline void SetEstimateVelocityX(int n, real v) {
    Vx[n]=v;
  }
  inline void SetProb(int n, real p) {
    Prob[n]=p;
  }
  inline void SetVelocityY(real v) {
    Uy=v;
  }
  inline void SetVelocityX(real v) {
    Ux=v;
  }
  inline void SetX(real x) {
    X=x;
  }
  inline void SetY(real y) {
    Y=y;
  }

  inline void SetNumPoints(int n, real link) {
    NumPoints=n;
	delete [] Link;
	Link=new real[NumPoints];
 	for(int i=0; i<NumPoints; ++i)
	   Link[i]=link;
  }

  inline void SetNumCandidates(int n) {
	  if(NumCandidates==n)
		  return;
    NumCandidates=n;
	if(Prob)
		delete [] Prob;
	if(Vy)
		delete [] Vy;
	if(Vx)
		delete [] Vx;
	Prob=new real[NumCandidates];
	Vy=new real[NumCandidates];
	Vx=new real[NumCandidates];
  }

  inline void SetLink(int i, real link) {
    Link[i]=link;
  }

  inline void Move(real rate, real var) {
	Y+=rate*(Uy+var*(rndm(0)-.5));
	X+=rate*(Ux+var*(rndm(0)-.5));
  }

  // Artithmetic overloaded operators

  // Standard Output

 protected:
  real* Vx; // velocity estimate
  real* Vy;
  real* Prob;
  real* Link;
  real X;
  real Y;
  real Ux;  // true velocity
  real Uy;
  int NumPoints;
  int NumCandidates;
};

ostream& operator<<(ostream& s, const Mpoint& v) {
    s << v.GetX() << ", " << v.GetY();
    s << ", " << v.GetVelocityX() << ", " << v.GetVelocityY() << endl;
	int i;
	s << "Link Weight:" << endl;
	for(i=0; i<v.GetNumPoints(); ++i)
		s << v.GetLink(i) << ", ";
	s << endl;
	s << "Velocity:" << endl;
	for(i=0; i<v.GetNumCandidates(); ++i) {
		s << v.GetEstimateVelocityX(i) << ", ";
		s << v.GetEstimateVelocityY(i) << ", ";
		s << v.GetProb(i) << "\n";
	}
	s << endl;
    return s;
}
#endif /* Mpoint_h */
