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
	bool bError = false;
	for (int i = 0; !bError && i < points.size(); ++i)
	{
		PointLight* w = points[i]->nextWinner();
		if (w != NULL)
		{
			if (points[i] == w->prevWinner()) 
			{
				merge(nodes[i], nodes[pmap[w]]);
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
			if (p->frame == q->frame)
			{
				if (findset(nodes[i]) == findset(nodes[j]))
				{
					p->print(NULL, "\n");
					q->print(">> ", "\n");
				}
			}
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
			P.push_back(new PointLight(x, y, fr, gid));
		}
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

	for (int i = 0; i < P.size(); ++i)
	{
		P[i]->initialize(P);
	}

	for (int i = 0; i<numIter; ++i) {
		printf("Iteration %d: \n", i);
		for (int i = 0; i < P.size(); ++i)
		{
			P[i]->update0(sigma, P[i]->Threshold, rate);
		}
		for (int i = 0; i < P.size(); ++i)
		{
			P[i]->update();
		}
	}

	vector<int> labels = clusterPoints(P);

	if (nlhs >= 1)
	{
		const int dims[] = { P.size(), 4 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			PointLight* p = P[i];
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], p->frame);
			SetData2(F, i, 3, dims[0], dims[1], (float)p->id);
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
		int nc = P[0]->NumCandidates;
		const int dims[] = { P.size(), 4 + 2 * nc + 2 * nc };
		vector<float> F(dims[0] * dims[1], -1.0f);
		for (int i = 0; i < dims[0]; ++i)
		{
			PointLight* p = P[i];
			SetData2(F, i, 0, dims[0], dims[1], (float)p->id);
			SetData2(F, i, 1, dims[0], dims[1], p->x);
			SetData2(F, i, 2, dims[0], dims[1], p->y);
			SetData2(F, i, 3, dims[0], dims[1], p->frame);
			int k = 4;
			for (int j = 0; j < nc; ++j)
			{
				if (j < p->prevC.size())
				{
					SetData2(F, i, k++, dims[0], dims[1], (float)p->prevC[j]->id);
					SetData2(F, i, k++, dims[0], dims[1], (float)p->prevP[j]);
				}
			}
			k = 4 + 2 * nc;
			for (int j = 0; j < nc; ++j)
			{
				if (j < p->nextC.size())
				{
					SetData2(F, i, k++, dims[0], dims[1], (float)p->nextC[j]->id);
					SetData2(F, i, k++, dims[0], dims[1], (float)p->nextP[j]);
				}
			}
		}
		plhs[2] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}

	for (int i = 0; i < P.size(); ++i)
	{
		delete P[i];
	}
	mexUnlock();
}

