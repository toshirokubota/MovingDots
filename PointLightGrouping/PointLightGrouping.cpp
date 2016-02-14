/*
This routine implements grouping of moving points.
*/
#include <mex.h>

#include <iostream>
#include <fstream>
#include <set>
#include <map>
using namespace std;
#include <stdlib.h>
#include <MCpoint.h>
#include <mexFileIO.h>
#include <szmexutilitytemplate.h>
#include <szMexUtility.h>
#include <szMiscOperations.h>
#include <PointLight.h>
#include <DisjointSet.h>

int PointLight::_id = 0;

vector<int>
clusterPoints(vector<PointLight*>& points)
{
	vector<Node<PointLight*>*> nodes;
	map<PointLight*, int> pmap;
	for (int i = 0; i < points.size(); ++i)
	{
		nodes.push_back(makeset(points[i]));
		pmap[points[i]] = i;
	}
	for (int i = 0; i < points.size(); ++i)
	{
		PointLight::PointLightTriple w = points[i]->winner();
		if (w.p != NULL)
		{
			PointLight::PointLightTriple wp = w.p->winner();
			if (wp.r == points[i]) 
			{
				merge(nodes[i], nodes[pmap[w.p]]);
			}
		}
		if (w.r != NULL)
		{
			PointLight::PointLightTriple wr = w.r->winner();
			if (wr.p == points[i])
			{
				merge(nodes[i], nodes[pmap[w.r]]);
			}
		}
	}
	vector<Node<PointLight*>*> reps = clusters(nodes);
	map<PointLight*, int> cmap;
	for (int i = 0; i < reps.size(); ++i)
	{
		cmap[reps[i]->key] = i;
	}
	vector<int> labels(points.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = cmap[findset(nodes[i])->key];
		labels[i] = k;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		PointLight* p = nodes[i]->key;
		for (int j = i + 1; j < nodes.size(); ++j)
		{
			PointLight* q = nodes[j]->key;
		}
		int k = cmap[findset(nodes[i])->key];
		labels[i] = k;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return labels;
}

vector<vector<PointLight*>>
organizeIntoFrames(vector<PointLight*>& P)
{
	vector<pair<float, PointLight*>> pairs(P.size());
	for (int i = 0; i < P.size(); ++i)
	{
		pairs[i].first = P[i]->frame;
		pairs[i].second = P[i];
	}
	sort(pairs.begin(), pairs.end());
	vector<vector<PointLight*>> frames;
	vector<PointLight*> frame;
	for (int i = 0; i < pairs.size(); ++i)
	{
		if (frame.empty() || frame[0]->frame == pairs[i].first)
		{
			frame.push_back(pairs[i].second);
		}
		else
		{
			frames.push_back(frame);
			frame.clear();
			frame.push_back(pairs[i].second);
		}
	}
	frames.push_back(frame);
	return frames;
}

/*
Check if there is already a point that is closely aligned to the line p-q and is temporarily closer to p.
*/
bool
isRedundant(PointLight* p, PointLight* q, vector<PointLight*>& Q, float thres)
{
	CParticleF p0(p->x, p->y);
	CParticleF q0(q->x, q->y);
	for (int i = 0; i < Q.size(); ++i)
	{
		CParticleF x(Q[i]->x, Q[i]->y);
		CParticleF c = Closest2Line(p0, q0, x);
		if (Distance(x, c) < thres)
		{
			return true;
		}
	}
	return false;
}

void
initialize(vector<PointLight*>& points, int range, float radius)
{
	vector<vector<PointLight*>> frames = organizeIntoFrames(points);
	float thres = radius / 10.0f;
	for (int i = 0; i < frames.size(); ++i)
	{
		for (int j = 0; j < frames[i].size(); j++)
		{
			vector<PointLight*> prev;
			PointLight* p = frames[i][j];
			for (int k = i - 1; k >= Max(0, i - range); k--)
			{
				vector<PointLight*> Q;
				for (int m = 0; m < frames[k].size(); ++m)
				{
					PointLight* q = frames[k][m];
					if (Distance(p->x, p->y, q->x, q->y) < radius)
					{
						if (isRedundant(p, q, prev, thres) == false)
						{
							Q.push_back(q);
						}
					}
				}
				prev.insert(prev.end(), Q.begin(), Q.end());
				if (prev.size() >= p->MaxNumNeighbors) break;
			}
			if (prev.empty()) prev.push_back(p);

			vector<PointLight*> next;
			for (int k = i + 1; k <= Min(frames.size() - 1, i + range); k++)
			{
				vector<PointLight*> Q;
				for (int m = 0; m < frames[k].size(); ++m)
				{
					PointLight* q = frames[k][m];
					if (Distance(p->x, p->y, q->x, q->y) < radius)
					{
						if (isRedundant(p, q, next, thres) == false)
						{
							Q.push_back(q);
						}
					}
				}
				next.insert(next.end(), Q.begin(), Q.end());
				if (next.size() >= p->MaxNumNeighbors) break;
			}
			if (next.empty()) next.push_back(p);

			for (int j = 0; j < prev.size(); ++j)
			{
				for (int k = 0; k < next.size(); ++k)
				{
					p->candidates.push_back(PointLight::Candidate(prev[j], p, next[k]));
				}
			}
			p->initializeProb();
		}
	}
}

void
update(vector<PointLight*>& P)
{
	for (int i = 0; i < P.size(); ++i)
	{
		if (P[i]->candidates.empty() == false)
		{
			//P[i]->updateParams();
			P[i]->updateFitness();
			//P[i]->updateProb();
		}
	}
	for (int i = 0; i < P.size(); ++i)
	{
		if (P[i]->candidates.empty() == false)
		{
			P[i]->updateProb();
		}
	}
}

vector<PointLight*>
removeDuplicates(vector<PointLight*>& P, float thres)
{
	vector<PointLight*> keep;
	for (int i = P.size() - 1; i >= 0; i--)
	{
		PointLight* p = P[i];
		bool bKeep = true;
		for (int j = i - 1; j >= 0; j--)
		{
			PointLight* q = P[j];
			if (p->frame == q->frame && Distance(p->x, p->y, q->x, q->y) < thres)
			{
				printf("remove %d(%f, %f)\n", p->id, p->x, p->y);
				bKeep = false;
				break;
			}
		}
		if (bKeep)
		{
			keep.insert(keep.begin(), p);
		}
	}
	return keep;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	printf("%s: This build was compiled at %s %s\n", "GroupPoints", __DATE__, __TIME__);
	if (nrhs < 0 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [A B C D] = GroupPoints(numiter, noise_var)");
		return;
	}
	PointLight::_id = 0;

	vector<PointLight*> P;
	{
		vector<float> P0;
		const int* dimsP;
		mxClassID classIdP;
		int ndimP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		for (int i = 0; i < dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			float fr = GetData2(P0, i, 2, dimsP[0], dimsP[1], (float)0);
			float gid = GetData2(P0, i, 3, dimsP[0], dimsP[1], (float)0);
			PointLight* pl = new PointLight(x, y, fr, gid);
			P.push_back(pl);
		}
		//P = removeDuplicates(P, 5.0);
	}

	int numIter = 1;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(numIter, prhs[1], classMode);
	}
	float sigma = 10.0;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(sigma, prhs[2], classMode);
	}
	float rate = 0.25;

	initialize(P, 5, 70.0);
	for (int i = 0; i < numIter; ++i)
	{
		printf("Iteration %d\n", i + 1);
		update(P);
	}
	vector<int> labels = clusterPoints(P);

	if (nlhs >= 1)
	{
		const int dims[] = { P.size(), 5 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			PointLight* p = P[i];
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], p->frame);
			SetData2(F, i, 3, dims[0], dims[1], (float)p->id);
			SetData2(F, i, 4, dims[0], dims[1], (float)p->gid);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		const int dims[] = { labels.size(), 1 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], labels[i]);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		int nc = P[0]->MaxNumCandidates;
		const int dims[] = { P.size(), nc};
		vector<int> F(dims[0] * dims[1], -1);
		for (int i = 0; i < dims[0]; ++i)
		{
			PointLight* p = P[i];
			for (int j = 0, k=0; j < P[i]->candidates.size(); ++j)
			{
				SetData2(F, i, k++, dims[0], dims[1], P[i]->candidates[j].tp.p->id);
				SetData2(F, i, k++, dims[0], dims[1], P[i]->candidates[j].tp.q->id);
				SetData2(F, i, k++, dims[0], dims[1], P[i]->candidates[j].tp.r->id);
			}
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}

	for (int i = 0; i < P.size(); ++i)
	{
		delete P[i];
	}
	mexUnlock();
}

