/*
This routine implements grouping of moving points.
*/
#include <mex.h>

#include <iostream>
#include <fstream>
#include <set>
using namespace std;
#include <stdlib.h>
#include <MCpoint.h>
#include <mexFileIO.h>
#include <szmexutilitytemplate.h>
#include <szMexUtility.h>
#include <szMiscOperations.h>
#include <DisjointSet.h>

int Factor=4;
int Width=256;
int Height=256;
float FrameRate=1.0;
float Rate=0.25;
//int NumIter=20;
//float NoiseVar=1.0;
int NumGroup=3;
int PointsPerGroup=8;
int NoisePoints=20;
int NumPoints=NumGroup*PointsPerGroup+NoisePoints;
float Epsilon=0.00001;
float Threshold1=.01;
float Threshold2=.01;
int Seed=54;
//float Sigma1=sqrt(2.0)*NoiseVar;
//float Sigma2=1.0;
//float Sigma3=1.0;

const int MaxNumGroups=10;
//const int MaxCandidates=10;
float MaxRange=25.0;
const int MaxNumPoints = 256;

float VelocityY[]={20.0, 20.0, -20.0};
float VelocityX[]={-20.0, 20.0, 20.0};

char infile[256];
char outfile[256];



vMCpoint
MovePoints(const vMCpoint& points, float rate, float var) 
{
	int nump=points.size();
	vMCpoint res=points;

	int i;
	for(i=0; i<nump; ++i) 
	{ 
		res[i].Move(rate,var);
	}
	return res;
}

vMCpoint
ShufflePoints(const vMCpoint& points)
{
	int nump=points.size();
	vMCpoint res(nump);

	int i;
	for(i=0; i<nump; ++i)
	{ 
		res[nump-i-1]=points[i];
	}
	return res;
}

float
Dist(float y1, float x1, float y2, float x2)
{
	return (y1-y2)*(y1-y2)+(x1-x2)*(x1-x2);
}

void
PlacePoints(vMCpoint& points, int height, int width) 
{
	int nump=points.size();
	int i,j,k,m;
	int offy=30;
	int offx=30;
	for(i=0,m=0; i<NumGroup; ++i)
	{
		for(j=0; j<PointsPerGroup; ++j,m++)
		{
			float y;
			float x;
			points[m].SetNumPoints(nump);
			y=(height-2*offy)*rndm(0)+offy;
			x=(width-2*offx)*rndm(0)+offx;
			points[m].SetY(y);
			points[m].SetX(x);
			points[m].SetVelocityY(VelocityY[i]);
			points[m].SetVelocityX(VelocityX[i]);
		}
	}
	for(i=NumGroup*PointsPerGroup; i<NumPoints; ++i)
	{
		points[i].SetNumPoints(nump);
		points[i].SetY(rndm(0)*height);
		points[i].SetX(rndm(0)*width);
		points[i].SetVelocityY(rndm(0)*50);
		points[i].SetVelocityX(rndm(0)*50);
	}
}

/*
Sort y, x and index according to the dist using bubble sort.
*/
void
SortCandidates(float* dist, float* y, float* x, int* index, int count) 
{
	for(int i=0; i<count; ++i) 
	{
		float max=dist[i];
		int maxi=i;
		for(int j=i; j<count; ++j)
		{
			if(max<dist[j])
			{
				maxi=j;
				max=dist[i];
			}
		}
		if(maxi!=i) 
		{
			float temp=dist[i];
			dist[i]=dist[maxi];
			dist[maxi]=temp;
			temp=y[i];
			y[i]=y[maxi];
			y[maxi]=temp;
			temp=x[i];
			x[i]=x[maxi];
			x[maxi]=temp;
			int tempi=index[i];
			index[i]=index[maxi];
			index[maxi]=tempi;
		}
	}
}

void
EstimateTransform(vMCpoint& p1, const vMCpoint& p2) 
{
	int nump=p1.size();
	float dist[MaxNumPoints];
	float est_vy[MaxNumPoints];
	float est_vx[MaxNumPoints];
	int index[MaxNumPoints];
	int i,j,k;
	for(i=0; i<nump; ++i) 
	{
		float x0=p1[i].GetX();
		float y0=p1[i].GetY();
		int count=0;
		for(j=0; j<nump; ++j) 
		{
			float x1=p2[j].GetX();
			float y1=p2[j].GetY();
			float dd=(x0-x1)*(x0-x1)+(y0-y1)*(y0-y1);
			if(dd<2.*MaxRange*MaxRange)
			{
				est_vy[count]=y1-y0;
				est_vx[count]=x1-x0;
				index[count]=j;
				count++;
			}
		}
		if(count==0) 
			cout << "Count=0 @ " << i << endl;
		cout << i << ":" << count << endl;
		//assert(count>0);
		SortCandidates(dist,est_vy,est_vx,index,count);
		count=Min(count,MCpointMaxCandidates);
		p1[i].SetNumCandidates(count);
		for(k=0; k<count; ++k) 
		{
			p1[i].SetEstimateVelocityY(k,est_vy[k]);
			p1[i].SetEstimateVelocityX(k,est_vx[k]);
			p1[i].SetIndex(k,index[k]);
			p1[i].SetProb(k,1.0/(count+1));
		}
		p1[i].SetProb(count, 1.0 / (count + 1));
	}
	//cout << "Done estimate transform.\n";
}

float
CompatibilityTransform(const MCpoint& p1, const MCpoint& p2, int m, int n, float sgm)
{
	float vy1=p1.GetEstimateVelocityY(m);
	float vx1=p1.GetEstimateVelocityX(m);
	float vy2=p2.GetEstimateVelocityY(n);
	float vx2=p2.GetEstimateVelocityX(n);
	float dd=Dist(vy1,vx1,vy2,vx2);
	if (dd > 2.*sgm*sgm)
	{
		return .0;
	}
	else
	{
		return exp(-dd / (2.*sgm*sgm));
	}
}

float
InterframeCompatibilityTransform(const MCpoint& p1, const MCpoint& p2, int m, int n, float sgm)
{
	float vy1=p1.GetEstimateVelocityY(m);
	float vx1=p1.GetEstimateVelocityX(m);
	float vy2=-p2.GetEstimateVelocityY(n);
	float vx2=-p2.GetEstimateVelocityX(n);
	float dd=Dist(vy1,vx1,vy2,vx2);
	if (dd > 2.*sgm*sgm)
	{
		return .0;
	}
	else
	{
		return exp(-dd / (2.*sgm*sgm));
	}
}

vMCpoint
UpdateTransformEstimate(const vMCpoint& points, float sigma, float rate)
{
	int nump=points.size();
	vMCpoint res=points;
	float sum[MCpointMaxCandidates];
	int i,m,n;
	//float nvy=0;
	//float nvx=0;
	float alpha=0.5;
	for(i=0; i<nump; ++i) 
	{
		float total_sum=0;
		MCpoint pnt1=points[i];
		for(m=0; m<pnt1.GetNumCandidates(); ++m) 
		{
			sum[m]=.0;
			float evy1=pnt1.GetEstimateVelocityY(m);
			float evx1=pnt1.GetEstimateVelocityX(m);
			float nvx = 0.0f;
			float nvy = 0.0f;
			for(int j=0; j<nump; ++j) 
			{
				if(i!=j) 
				{
					MCpoint pnt2=points[j];
					for(n=0; n<pnt2.GetNumCandidates(); ++n)
					{
						float evy2=pnt2.GetEstimateVelocityY(n);
						float evx2=pnt2.GetEstimateVelocityX(n);
						float lnk1=pnt1.GetLink(j,m,n);
						float lnk2=pnt2.GetLink(i,n,m);
						float dd=CompatibilityTransform(pnt1,pnt2,m,n,sigma);
						float dy=evy2-evy1;
						float dx=evx2-evx1;
						float pp=pnt2.GetProb(n);
						nvy+=dd*pp*lnk1*lnk2*dy;
						nvx+=dd*pp*lnk1*lnk2*dx;
						/*if (i == 9 && m == 2)
						{
							printf("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
								evx2, evy2, lnk1, lnk2, dd, dx, dy, pp, nvx, nvy);
						}*/
					}
				}
			}
			/*if (i == 9 && m == 2)
			{
				printf("Final: %f, %f, %f, %f, %f, %f\n",
					evx1, evy1, nvx, nvy, evx1 + rate*nvx, evy1 + rate*nvy);
			}*/

			res[i].SetEstimateVelocityY(m,evy1+rate*nvy);
			res[i].SetEstimateVelocityX(m,evx1+rate*nvx);
		}
	}
	return res;
}

vMCpoint
UpdateLinkWeights(const vMCpoint& points, float thres, float sigma) 
{
	int nump=points.size();
	vMCpoint res=points;
	for(int i=0; i<nump; ++i)
	{
		MCpoint p1=points[i];
		int numc1=p1.GetNumCandidates();
		for(int m=0; m<numc1; ++m) 
		{
			for(int j=0; j<nump; ++j) 
			{
				MCpoint p2=points[j];
				int numc2=p2.GetNumCandidates();
				for(int n=0; n<numc2; ++n) 
				{
					if(i!=j)
					{
						float dd=CompatibilityTransform(p1,p2,m,n,sigma);
						float lnk1=p1.GetLink(j,m,n);
						float lnk2=p2.GetLink(i,n,m);
						float pp=p2.GetProb(n);
						float ee=pp*lnk1*lnk2;
						float pp2=ee*dd/(ee*dd+(1.-ee)*thres);
						if (pp2 == pp2)
						{
							res[i].SetLink(j, m, n, pp2);
						}
						else
						{
							res[i].SetLink(j, m, n, 0.5f);
						}
					}
				}
			}
		}
	}
	return res;
}

vMCpoint
UpdateProbMeasure(const vMCpoint& frame1, const vMCpoint& frame2, 
				  float thres, float sigma)
{
	int nump=frame1.size();
	vMCpoint res=frame1;
	float sum[MCpointMaxCandidates];
	for(int i=0; i<nump; ++i)
	{
		MCpoint p1=frame1[i];
		int numc1=p1.GetNumCandidates();
		float total_sum=0;
		int m;
		for(m=0; m<numc1; ++m) 
		{ //compute fitness for each candidate
			sum[m]=.0;
			MCpoint p3=frame2[p1.GetIndex(m)];
			float match_prob=.0;
			//check if the candidate also consider it as a candidate
			for(int m2=0; m2<p3.GetNumCandidates(); ++m2) 
			{
				if(p3.GetIndex(m2)==i)
				{
					match_prob=p3.GetProb(m2);
					break;
				}
			}
			//if(match_prob>0) 
			{ //mutual candidancy
				for(int j=0; j<nump; ++j)
				{
					//compute the support from other INTRA-frame points
					//measured by 1. transformation similarity
					//            2. strength of links
					if(i!=j)
					{
						MCpoint p2=frame1[j];
						int numc2=p2.GetNumCandidates();
						for(int n=0; n<numc2; ++n)
						{
							float dd=CompatibilityTransform(p1,p2,m,n,sigma);
							float lnk1=p1.GetLink(j,m,n);
							float lnk2=p2.GetLink(i,n,m);
							float pp=p2.GetProb(n);
							float ss=pp*lnk1*lnk2;
							sum[m]+=ss*match_prob;
						}
					}
				}
				sum[m]*=p1.GetProb(m);
				total_sum+=sum[m];
			}
		}
		if (total_sum > 1.0e-5)
		{
			for (m = 0; m < numc1; ++m)
			{
				res[i].SetProb(m, sum[m] / (total_sum + p1.GetProb(numc1)*thres));
			}
			res[i].SetProb(numc1, p1.GetProb(numc1)*thres / (total_sum + p1.GetProb(numc1)*thres));
		}
		else
		{
			for (m = 0; m < numc1; ++m)
			{
				res[i].SetProb(m, 0.0f);
			}
			res[i].SetProb(numc1, 1.0f);
		}
	}
	return res;
}

vector<int>
labelIntraFramePoints(vector<MCpoint>& points)
{
	vector<Node<int>*> vnodes;
	for (int i = 0; i < points.size(); ++i)
	{
		vnodes.push_back(makeset(i));
	}
	for (int i = 0; i < points.size(); ++i)
	{
		for (int j = i + 1; j < points.size(); ++j)
		{
			//if (i == j) continue;
			for (int m = 0; m < points[i].GetNumCandidates(); ++m)
			{
				for (int n = 0; n<points[j].GetNumCandidates(); ++n)
				{
					if (points[i].GetLink(j, m, n) > .5) {
						merge(vnodes[i], vnodes[j]);
					}
				}
			}
		}
	}
	vector<Node<int>*> labels = clusters(vnodes);
	vector<int> ilabels(points.size());
	for (int i = 0; i < points.size(); ++i)
	{
		ilabels[i] = distance(labels.begin(), find(labels.begin(), labels.end(), findset(vnodes[i])));
	}
	for (int i = 0; i < points.size(); ++i)
	{
		delete vnodes[i];
	}
	return ilabels;
}

vector<pair<int, int>>
findIntraframeMatches(vector<MCpoint>& frame)
{
	vector<pair<int, int>> matches;
	for (int i = 0; i < frame.size(); ++i)
	{
		for (int j = 0; j < frame.size(); ++j)
		{
			for (int m = 0; m < frame[i].GetNumCandidates(); ++m)
			{
				for (int n = 0; n < frame[j].GetNumCandidates(); ++n)
				{
					float w = frame[i].GetLink(j, m, n);
					if (w > 0.5)
					{
						matches.push_back(pair<int, int>(i, j));
						break;
					}
				}
			}
		}
	}
	return matches;
}

vector<pair<int, int>>
findInterframeMatches(vector<MCpoint>& frame1, vector<MCpoint>& frame2)
{
	vector<pair<int, int>> matches;
	for (int i = 0; i < frame1.size(); ++i)
	{
		for (int m = 0; m < frame1[i].GetNumCandidates(); ++m)
		{
			float p = frame1[i].GetProb(m);
			if (p > 0.5)
			{
				int j = frame1[i].GetIndex(m);
				for (int n = 0; n < frame2[j].GetNumCandidates(); ++n)
				{
					float q = frame2[j].GetProb(n);
					if (q > 0.5 && frame2[j].GetIndex(n) == i)
					{
						matches.push_back(pair<int, int>(i, j));
						break;
					}
				}
			}
		}
	}
	return matches;
}

bool sanityCheck(vector<MCpoint>& points)
{
	for (int i = 0; i < points.size(); ++i)
	{
		for (int j = 0; j < points[i].GetNumCandidates(); ++j)
		{
			float p = points[i].GetProb(j);
			if (p != p)
			{
				printf("sanityCheck: prob = nan on [%d, %d]\n", i, j);
				return false;
			}
			if (p < 0 || p>1.0f)
			{
				printf("sanityCheck: prob = %f on [%d, %d]\n", p, i, j);
				return false;
			}
			float vx = points[i].GetEstimateVelocityX(j);
			float vy = points[i].GetEstimateVelocityY(j);
			if (vx != vx)
			{
				printf("sanityCheck: vx nan on [%d, %d]\n", i, j);
				return false;
			}
			if (vy != vy)
			{
				printf("sanityCheck: vy nan on [%d, %d]\n", i, j);
				return false;
			}
		}
	}
	return true;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "GroupPoints", __DATE__, __TIME__);
	if (nrhs < 0 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [A B C D] = GroupPoints(numiter, noise_var)");
		return;
	}
	int numIter = 20;
	if (nrhs >= 1)
	{
		mxClassID classMode;
		ReadScalar(numIter, prhs[0], classMode);
	}
	float noiseVar = 1.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(noiseVar, prhs[1], classMode);
	}
	float sigma = sqrt(2.0) * noiseVar;

	vMCpoint frame1(NumPoints);
	vMCpoint frame2;
	vMCpoint tmp_frm;
	PlacePoints(frame1,Height,Width);
	frame2 = MovePoints(frame1, FrameRate, noiseVar);
	frame2=ShufflePoints(frame2); //to make sure the array index is not helping
	EstimateTransform(frame1,frame2);
	EstimateTransform(frame2,frame1);

	for (int i = 0; i<numIter; ++i) {
		printf("Iteration %d: \n", i);
		frame1 = UpdateLinkWeights(frame1,Threshold1,sigma);
		frame2 = UpdateLinkWeights(frame2, Threshold1, sigma);
		frame1 = UpdateTransformEstimate(frame1, sigma, Rate);
		frame2 = UpdateTransformEstimate(frame2, sigma, Rate);
		tmp_frm = UpdateProbMeasure(frame1, frame2, Threshold2, sigma);
		frame2 = UpdateProbMeasure(frame2, tmp_frm, Threshold2, sigma);
		frame1=tmp_frm;
		if (sanityCheck(frame1) == false)
		{
			printf("frame1 - failed at sanity check at %d\n", i);
			//mexErrMsgTxt("Sanity check failed on frame1.\n");
			break;
		}
		if (sanityCheck(frame2) == false)
		{
			printf("frame2 - failed at sanity check at %d\n", i);
			//mexErrMsgTxt("Sanity check failed on frame2.\n");
			break;
		}
	}
	printf("Frame 1:\n");
	for (int i = 0; i<frame1.size(); ++i)
	{
		printf("Point #%d\n", i);
		frame1[i].print();
	}
	printf("Frame 2:\n");
	for (int i = 0; i<frame2.size(); ++i)
	{
		printf("Point #%d\n", i);
		frame2[i].print();
	}
	vector<int> vlabel1 = labelIntraFramePoints(frame1);
	vector<int> vlabel2 = labelIntraFramePoints(frame2);
	if (nlhs >= 1)
	{
		const int dims[] = { frame1.size(), 3 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], frame1[i].GetX());
			SetData2(F, i, 1, dims[0], dims[1], frame1[i].GetY());
			SetData2(F, i, 2, dims[0], dims[1], (float)vlabel1[i]);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);

	}
	if (nlhs >= 2)
	{
		const int dims[] = { frame2.size(), 3 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], frame2[i].GetX());
			SetData2(F, i, 1, dims[0], dims[1], frame2[i].GetY());
			SetData2(F, i, 2, dims[0], dims[1], (float)vlabel2[i]);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);

	}
	if (nlhs >= 3)
	{
		vector<pair<int, int>> matches = findIntraframeMatches(frame1);
		const int dims[] = { matches.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], matches[i].first+1);
			SetData2(F, i, 1, dims[0], dims[1], matches[i].second + 1);
			SetData2(F, i, 2, dims[0], dims[1], vlabel1[matches[i].first]);
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 4)
	{
		vector<pair<int, int>> matches = findIntraframeMatches(frame2);
		const int dims[] = { matches.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], matches[i].first+1);
			SetData2(F, i, 1, dims[0], dims[1], matches[i].second+1);
			SetData2(F, i, 2, dims[0], dims[1], vlabel2[matches[i].first]);
		}
		plhs[3] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 5)
	{
		vector<pair<int, int>> matches = findInterframeMatches(frame1, frame2);
		const int dims[] = { matches.size(), 4 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], matches[i].first+1);
			SetData2(F, i, 1, dims[0], dims[1], matches[i].second + 1);
			SetData2(F, i, 2, dims[0], dims[1], vlabel1[matches[i].first]);
			SetData2(F, i, 3, dims[0], dims[1], vlabel2[matches[i].second]);
		}
		plhs[4] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	mexUnlock();
}

