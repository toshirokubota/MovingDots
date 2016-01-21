#ifndef _MotionPoint_h_
#define _MotionPoint_h_

#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;
#include <szMexUtility.h>
#include <Feature.h>
#include <mex.h>

struct PointLight
{
	const float Threshold = 0.05;
	const int NumCandidates = 3;
	PointLight(float x0 = 0, float y0 = 0, float fr = 0)
	{
		x = x0;
		y = y0;
		frame = fr;
		id = _id++;
	}

	void initialize(vector<PointLight*>& points);
	void update0(float sigma, float thres, float rate);
	void update()
	{
		for (int i = 0; i < prevP.size(); ++i)
		{
			prevP[i] = prevP0[i];
		}
		for (int i = 0; i < nextP.size(); ++i)
		{
			nextP[i] = nextP0[i];
		}
	}
	void print(char* tab=NULL, char* newl=NULL);

	float x;
	float y;
	float frame; 
	int id;

	vector<PointLight*> prevC; //possible correspondences from prev frame
	vector<PointLight*> nextC; //possible correspondences from prev frame
	vector<float> prevP;
	vector<float> nextP;

	vector<float> prevP0; //temporary copy of inter-framce correspondence probability
	vector<float> nextP0; //temporary copy of inter-framce correspondence probability

	static int _id;
};

#endif /* _MotionPoint_H */
