/*
This routine implements grouping of moving points.
*/
#include <mex.h>

#include <iostream>
#include <fstream>
#include <set>
using namespace std;
#include <stdlib.h>
#include <MotionPoint.h>
#include <MotionPointProcessing.h>
#include <mexFileIO.h>
#include <szmexutilitytemplate.h>
#include <szMexUtility.h>
#include <szMiscOperations.h>
#include <DisjointSet.h>

int MotionPoint::_id = 0;
//int Factor = 4;
//int Width=256;
//int Height=256;
//float FrameRate=1.0;
float Rate=0.25;
//int NumIter=20;
//float NoiseVar=1.0;
//int NumGroup=3;
//int PointsPerGroup=8;
//int NoisePoints=20;
//int NumPoints=NumGroup*PointsPerGroup+NoisePoints;
//float Epsilon=0.00001;
float Threshold1=.01;
float Threshold2=.01;
//int Seed=54;
//float Sigma1=sqrt(2.0)*NoiseVar;
//float Sigma2=1.0;
//float Sigma3=1.0;

//const int MaxNumGroups=10;
//const int MaxCandidates=10;
float MaxRange=25.0;
//const int MaxNumPoints = 256;
const int NumMaxCandidates = 10;

//float VelocityY[]={20.0, 20.0, -0.0};
//float VelocityX[]={-20.0, 20.0, 20.0};

//char infile[256];
//char outfile[256];


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "GroupPoints", __DATE__, __TIME__);
	MotionPoint::_id = 0;

	if (nrhs < 2 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [A B C D] = MovingSkeletonGrouping(P, Q)");
		return;
	}
	//Points
	vector<MotionPoint*> frame1;
	{
		vector<float> P0;
		mxClassID classIdP;
		const int* dimsP;
		int ndimP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		for (int i = 0; i < dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			frame1.push_back(new MotionPoint(x, y));
		}
	}
	vector<MotionPoint*> frame2;
	{
		vector<float> P0;
		mxClassID classIdP;
		const int* dimsP;
		int ndimP;
		LoadData(P0, prhs[1], classIdP, ndimP, &dimsP);
		for (int i = 0; i < dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			frame2.push_back(new MotionPoint(x, y));
		}
	}
	int numIter = 10;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(numIter, prhs[2], classMode);
	}
	int numNoise = 20;
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(numNoise, prhs[3], classMode);
	}
	float sigma = 10.0f;
	if (nrhs >= 5)
	{
		mxClassID classMode;
		ReadScalar(sigma, prhs[4], classMode);
	}
	//float sigma = Max(1.0f, sqrt(2.0) * noiseVar);

	EstimateTransform(frame1, frame2, NumMaxCandidates, MaxRange);
	EstimateTransform(frame2, frame1, NumMaxCandidates, MaxRange);

	for (int i = 0; i < frame1.size(); ++i)
	{
		frame1[i]->initializeLinks(frame1);
	}
	for (int i = 0; i < frame2.size(); ++i)
	{
		frame2[i]->initializeLinks(frame2);
	}

	for (int i = 0; i<numIter; ++i) {
		printf("Iteration %d: \n", i);
		UpdateLinkWeights(frame1,Threshold1,sigma);
		UpdateLinkWeights(frame2, Threshold1, sigma);
		UpdateTransformEstimate(frame1, sigma, Rate);
		UpdateTransformEstimate(frame2, sigma, Rate);
		UpdateProbMeasureTmp(frame1, frame2, Threshold2, sigma);
		UpdateProbMeasureTmp(frame2, frame1, Threshold2, sigma);
		UpdateProbMeasure(frame1);
		UpdateProbMeasure(frame2);
	}
	vector<int> vlabel1 = labelIntraFramePoints(frame1);
	vector<int> vlabel2 = labelIntraFramePoints(frame2);
	if (nlhs >= 1)
	{
		const int dims[] = { frame1.size(), 6 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			MotionPointSynth* p = (MotionPointSynth*) frame1[i];
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], (float)vlabel1[i]);
			SetData2(F, i, 3, dims[0], dims[1], p->ux);
			SetData2(F, i, 4, dims[0], dims[1], p->uy);
			SetData2(F, i, 5, dims[0], dims[1], (float)p->groupId);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);

	}
	if (nlhs >= 2)
	{
		const int dims[] = { frame2.size(), 6 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			MotionPointSynth* p = (MotionPointSynth*)frame2[i];
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], (float)vlabel2[i]);
			SetData2(F, i, 3, dims[0], dims[1], p->ux);
			SetData2(F, i, 4, dims[0], dims[1], p->uy);
			SetData2(F, i, 5, dims[0], dims[1], (float)p->groupId);
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
			SetData2(F, i, 0, dims[0], dims[1], matches[i].first + 1);
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

	for (int i = 0; i < frame1.size(); ++i)
	{
		delete frame1[i];
	}
	for (int i = 0; i < frame2.size(); ++i)
	{
		delete frame2[i];
	}
	mexUnlock();
}

