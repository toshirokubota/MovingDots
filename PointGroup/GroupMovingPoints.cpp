/*
This routine implements grouping of moving points.
*/

#include <iostream.h>
//#include <fstream.h>
//using namespace std;

#include <mytype.h>
#include <ImageProc.h>
#include <Feature.h>
#include <Convolve.h>
#include <Gradient.h>
#include <Misc.h>

#include <Draw.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <Mpoint.h>

int Factor=4;
int Width=256;
int Height=256;
real FrameRate=1.0;
real Rate=0.25;
int NumIter=20;
real NoiseVar=1.0;
const int NumGroup=3;
const int PointsPerGroup=3;
int NumPoints=NumGroup*PointsPerGroup;
real Epsilon=0.00001;
real Threshold=.1;
int Seed=0;
real Sigma1=1.0;
real Sigma2=1.0;

const int NumGroups=3;
const int MaxCandidates=10;

real VelocityY[NumGroup]={0.0, 2.0, -2.0};
real VelocityX[NumGroup]={0.0, 2.0, 2.0};

char infile[256];
char outfile[256];

//
// Usage - argument parsing routines
//
void
Usage(char* program) {
	cerr << program << ": usage" << endl;
	cerr << "Required arguments:" << endl;
	//cerr << "\t -i <input file name>" << endl;
	cerr << "\t -o <output file name>" << endl;
	
	cerr << "Optional arguments:" << endl;
	cerr << "\t -m (int) set the # of iterations" << endl;
	cerr << "\t \t default value is " << NumIter << endl;
	cerr << "\t -rt (int) set the update rate (>0)" << endl;
	cerr << "\t \t default value is " << Rate << endl;
	cerr << "\t -nv (int) set the noise variance" << endl;
	cerr << "\t \t default value is " << NoiseVar << endl;
	cerr << "\t -np (int) set the # of sample points" << endl;
	cerr << "\t \t default value is " << NumPoints << endl;
	cerr << "\t -seed (int) set the seed value for random number generator" << endl;
	cerr << "\t \t default value is " << Seed << endl;
}

void
PrintParameters() {
	cout << "input file: " << infile << endl;
	cout << "output file: " << outfile << endl;
	cout << "# of iterations: " << NumIter << endl;
	cout << "the update rate: " << Rate << endl;
	cout << "the noise variance: " << NoiseVar << endl;
	cout << "# of sample points: " << NumPoints << endl;
	cout << "Seed number: " << Seed << endl;
}


int
ReadArguments(int argc, char* argv[]) {
	int in_found=0;
	int out_found=0;
	real tempval;
	for(int i=1; i<argc; ++i) {
		if(strcmp(argv[i],"-i")==0) {
			if(argc>i+1) {
				in_found=1;
				strcpy(infile,argv[++i]);
			}
			else
				return 1; // no file specified
		}
		else if(strcmp(argv[i],"-o")==0) {
			if(argc>i+1) {
				out_found=1;
				strcpy(outfile,argv[++i]);
			}
			else
				return 2; // no file specified
		}
		else if(strcmp(argv[i],"-m")==0) {
			if(argc>i+1) {
				NumIter=atoi(argv[++i]);
			}
			else
				return 6; // no value specified
		}
		else if(strcmp(argv[i],"-seed")==0) {
			if(argc>i+1) {
				Seed=atoi(argv[++i]);
			}
			else
				return 6; // no value specified
		}
		else if(strcmp(argv[i],"-np")==0) {
			if(argc>i+1) {
				NumPoints=atoi(argv[++i]);
			}
			else
				return 6; // no value specified
		}
		else if(strcmp(argv[i],"-rt")==0) {
			if(argc>i+1) {
				tempval=strtod(argv[++i],NULL);
				if(tempval<0)
					return 8;
				else
					Rate=tempval;
			}
			else
				return 6; // no value specified
		}
		else if(strcmp(argv[i],"-nv")==0) {
			if(argc>i+1) {
				NoiseVar=strtod(argv[++i],NULL);
			}
			else
				return 6; // no value specified
		}
	}
	
	//if(!in_found)
	//  return 11;
	rndm(Seed);
	if(!out_found)
		return 12;
	else
		return 0;
}

void
CopyPoints(Mpoint* dest, Mpoint* src, int nump) {
	for(int i=0; i<nump; ++i) 
		dest[i]=src[i];
}

RealImage
ShowPoints(Mpoint* points, int h, int w, int nump) {
	RealImage im(1,h,w,0);
	int i;
	int length=5;
	for(i=0; i<nump; ++i) {
		real x=points[i].GetX();
		real y=points[i].GetY();
		int ix=Round(x);
		int iy=Round(y);
		if(iy>=0 && iy<h && ix>=0 && ix<w)
			im.SetPixel(0,iy,ix,1.0);
		//DrawCross(im,0,iy,ix,0);
	}
	return im;
}

RealImage
ShowLinks(Mpoint* points, int h, int w, int nump) {
	RealImage im(1,h,w,0);
	int i,j;
	int length=5;
	for(i=0; i<nump; ++i) {
		real x=points[i].GetX();
		real y=points[i].GetY();
		int ix=Round(x);
		int iy=Round(y);
		for(j=0; j<nump; ++j) {
			if(i!=j) {
				real x2=points[j].GetX();
				real y2=points[j].GetY();
				int ix2=Round(x2);
				int iy2=Round(y2);
				real ww=points[i].GetLink(j);
				if(ww>0.1)
					DrawLine(im,0,iy,ix,iy2,ix2,ww);
			}
		}
	}
	return im;
}

void
MovePoints(Mpoint* points, real rate, int nump, real var) {
	
	int i;
	for(i=0; i<nump; ++i) {
		points[i].Move(rate,var);
	}
}

real
Dist(real y1, real x1, real y2, real x2) {
	return (y1-y2)*(y1-y2)+(x1-x2)*(x1-x2);
}

void
InitializePoints(Mpoint* points, Mpoint* oldpoints, 
				 int height, int width, int nump, real rate, real var) {
	RealImage map(3,height,width,0);
	int i,j,k;
	int offy=30;
	int offx=30;
	for(i=0; i<nump; ++i) {
		real y=(height-2*offy)*rndm(0)+offy;
		real x=(width-2*offx)*rndm(0)+offx;
		int gid=i/PointsPerGroup;
		points[i].SetNumPoints(nump,.5);
		points[i].SetY(y);
		points[i].SetX(x);
		points[i].SetVelocityY(VelocityY[gid]);
		points[i].SetVelocityX(VelocityX[gid]);
		map.SetPixel(gid,(int)y,(int)x,1);
	}
	map.WritePnmFile("temp2.ppm",IO_PPM,1);
	CopyPoints(oldpoints,points,nump);
	MovePoints(points,rate,nump,var);
	
	real* est_vy=new real[MaxCandidates];
	real* est_vx=new real[MaxCandidates];
	for(i=0; i<nump; ++i) {
		real x0=oldpoints[i].GetX();
		real y0=oldpoints[i].GetY();
		int count=0;
		for(j=0; j<nump; ++j) {
			real x1=points[j].GetX();
			real y1=points[j].GetY();
			real dd=(x0-x1)*(x0-x1)+(y0-y1)*(y0-y1);
			if(dd<16.0 && count<MaxCandidates) {
				est_vy[count]=y1-y0;
				est_vx[count]=x1-x0;
				count++;
			}
		}
		if(count==0) 
			cout << "Count=0 @ " << i << endl;
		//assert(count>0);
		points[i].SetNumCandidates(count);
		cout << i << ":" << count << endl;
		for(k=0; k<count; ++k) {
			points[i].SetEstimateVelocityY(k,est_vy[k]);
			points[i].SetEstimateVelocityX(k,est_vx[k]);
			points[i].SetProb(k,1.0/count);
		}
	}
}

real
CompatibilityVelocity(const Mpoint& p1, const Mpoint& p2, int m, int n, real sgm) {
	real vy1=p1.GetEstimateVelocityY(m);
	real vx1=p1.GetEstimateVelocityX(m);
	real vy2=p2.GetEstimateVelocityY(n);
	real vx2=p2.GetEstimateVelocityX(n);
	real dd=Dist(vy1,vx1,vy2,vx2);
	if(dd>sgm)
		return .0;
	else
		return exp(-dd/(2.*sgm*sgm));
}

real
CompatibilityPoints(const Mpoint& point1, const Mpoint& point2, real thres, real sgm) {
	real maxc=.0;
	for(int i=0; i<point1.GetNumCandidates(); ++i) {
		real vy=point1.GetEstimateVelocityY(i);
		real vx=point1.GetEstimateVelocityX(i);
		real p1=point1.GetProb(i);
		if(p1<thres)
			continue;
		for(int j=0; j<point2.GetNumCandidates(); ++j) {
			real vy2=point2.GetEstimateVelocityY(j);
			real vx2=point2.GetEstimateVelocityX(j);
			real dd=Dist(vy,vx,vy2,vx2);
			real p2=point2.GetProb(j);
			if(p2<thres)
				continue;
			real cc=exp(-dd/(2.*sgm*sgm));
			maxc=Max(cc,maxc);
		}
	}
	return maxc;
}


void
UpdateVelocityEstimate(Mpoint* points, Mpoint* prev, int nump, real sigma, real rate) {
	real* sum=new real[MaxCandidates];
	int i,m,n;
	real nvy=0;
	real nvx=0;
	real alpha=0.5;
	for(i=0; i<nump; ++i) {
		real total_sum=0;
		for(m=0; m<points[i].GetNumCandidates(); ++m) {
			sum[m]=.0;
			real evy=prev[i].GetEstimateVelocityY(m);
			real evx=prev[i].GetEstimateVelocityX(m);
			for(int j=0; j<nump; ++j) {
				if(i!=j) {
					real p1=prev[i].GetLink(j);
					real p2=prev[j].GetLink(i);
					for(n=0; n<points[j].GetNumCandidates(); ++n) {
						real cc=CompatibilityVelocity(prev[i],prev[j],m,n,sigma);
						sum[m]+=p2*cc;
						real dy=prev[j].GetEstimateVelocityY(n)-evy;
						real dx=prev[j].GetEstimateVelocityX(n)-evx;
						nvy+=cc*p2*dy;
						nvx+=cc*p2*dx;
					}
				}
			}
			real vy=points[i].GetY()-prev[i].GetY();
			real vx=points[i].GetX()-prev[i].GetX();
			points[i].SetEstimateVelocityY(n,alpha*evy+(1.-alpha)*vy+rate*nvy);
			points[i].SetEstimateVelocityX(n,alpha*evx+(1.-alpha)*vx+rate*nvx);
			sum[m]*=prev[i].GetProb(m);
			total_sum+=sum[m];
		}
		for(m=0; m<points[i].GetNumCandidates(); ++m) {
			points[i].SetProb(m,sum[m]/total_sum);
		}
	}
}

void
UpdateLinkWeights(Mpoint* points, Mpoint* prev, int nump, real thres, real sigma) {
	for(int i=0; i<nump; ++i) {
		for(int j=0; j<nump; ++j) {
			if(i!=j) {
				real cc=CompatibilityPoints(prev[i],prev[j],0.1,sigma);
				real pp=prev[i].GetLink(j);
				real p2=pp*cc/(pp*cc+(1.-pp)*thres);
				points[i].SetLink(j,p2);
			}
		}
	}
}

void main(int argc,char *argv[]) {
	
	if(ReadArguments(argc,argv)!=0) {
		Usage(argv[0]);
		abort();
	}
	else
		PrintParameters();
	
	Mpoint* points=new Mpoint[NumPoints];
	Mpoint* previous=new Mpoint[NumPoints];
	InitializePoints(points,previous,Height,Width,NumPoints,FrameRate,NoiseVar);
	//CopyPoints(orig,points,NumPoints);
	RealImage map=ShowLinks(points,Height,Width,NumPoints);
	InteractiveImageWrite(map,"temp.pgm",1);
	
	for(int i=0; i<NumIter; ++i) {
		CopyPoints(previous,points,NumPoints);
		MovePoints(points,FrameRate,NumPoints,NoiseVar);
		UpdateVelocityEstimate(points,previous,NumPoints,Sigma1,Rate);
		UpdateLinkWeights(points,previous,NumPoints,Threshold,Sigma2);
		map=ShowLinks(points,Height,Width,NumPoints);
		sprintf(outfile,"out%2.2d.pgm",i);
		map.WritePnmFile(outfile,IO_PGM,1);
	}
	for(i=0; i<NumPoints; ++i)
		cout << points[i];
}