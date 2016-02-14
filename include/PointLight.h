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
			if (r != NULL && p != NULL && (p != q || r != q))
			{
				float t = r->frame - q->frame;
				float s = p->frame - q->frame;
				if (t == 0 && s == 0)
				{
					v[1] = v[2] = v[4] = v[5] = 0;
				}
				else if (t == 0)
				{
					v[1] = (p->x - q->x) / s;
					v[4] = (p->y - q->y) / s;
					v[2] = v[5] = 0;
				}
				else if (s == 0)
				{
					v[1] = (r->x - q->x) / t;
					v[4] = (r->y - q->y) / t;
					v[2] = v[5] = 0;
				}
				else
				{
					float denom = s * t * (s - t);
					v[1] = (q->x * (t*t - s*s) + s*s*r->x - t*t*p->x) / denom;
					v[2] = (q->x*(s - t) - s*r->x + t*p->x) / denom;
					v[4] = (q->y * (t*t - s*s) + s*s*r->y - t*t*p->y) / denom;
					v[5] = (q->y*(s - t) - s*r->y + t*p->y) / denom;
				}
			}
			else
			{
				v[1] = 0;
				v[2] = 0;
				v[4] = 0;
				v[5] = 0;
			}
		}
		CParticleF evaluate(float frame)
		{
			float t = frame - q->frame;
			float x = v[0] + t * v[1] + t*t*v[2];
			float y = v[3] + t * v[4] + t*t*v[5];
			return CParticleF(x, y);
		}
		CParticleF velocity(float frame)
		{
			float t = frame - q->frame;
			float x = v[1] + 2*t*v[2];
			float y = v[4] + 2*t*v[5];
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
	void updateParams();
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
	static void _FindOptimumParameters(Candidate& a, Candidate& b, float res[6], float G[3][3], float H[3][3]);
	static int _id;
};

#endif /* _MotionPoint_H */
