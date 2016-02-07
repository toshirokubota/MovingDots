#ifndef _MotionPoint_h_
#define _MotionPoint_h_

#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;
#include <szMexUtility.h>
#include <Feature.h>
#include <mex.h>
#include <szParticleF.h>

struct PointLight
{
	struct PointLightTriple
	{
		PointLightTriple(PointLight* p, PointLight* q, PointLight* r)
		{
			this->p = p;
			this->q = q;
			this->r = r;
			v[0] = q->x;
			v[3] = q->y;
			if (r != NULL && p != NULL)
			{
				v[1] = (r->x - p->x) / (r->frame - p->frame);
				v[4] = (r->y - p->y) / (r->frame - p->frame);
			}
			else
			{
				v[1] = 0;
				v[4] = 0;
			}
			v[2] = 0;
			v[5] = 0;
		}
		CParticleF evaluate(float frame)
		{
			float t = frame - q->frame;
			float x = v[0] + t * v[1];
			float y = v[3] + t * v[4];
			return CParticleF(x, y);
		}
		PointLight* p;
		PointLight* q;
		PointLight* r;
		float v[6];
	};
	struct Candidate
	{
		Candidate(PointLight* p, PointLight* q, PointLight* r) : tp(p, q, r)
		{
			comp = 0; // _Compatibility(p, q, sigma);
			prob = 0;
		}
		PointLightTriple tp;
		float comp; //compatibility
		float prob;
	};

	const float Threshold = 0.01;
	const int MaxNumNeighbors = 5;
	const int MaxNumCandidates = MaxNumNeighbors * MaxNumNeighbors;
	PointLight(float x0 = 0, float y0 = 0, float fr = 0, int g = 0)
	{
		x = x0;
		y = y0;
		frame = fr;
		id = _id++;
		gid = g;
	}

	void print(char* tab=NULL, char* newl=NULL);
	void updateFitness();
	void updateProb();
	void initializeProb();
	PointLightTriple winner();

	float x;
	float y;
	float frame; 
	int id;
	int gid; //group id for debuggin.

	vector<Candidate> candidates; //possible correspondences from prev frame

	static float _Compatibility(Candidate& a, Candidate& b);
	static int _id;
};

#endif /* _MotionPoint_H */
