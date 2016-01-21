#include <iostream>
using namespace std;
#include <stdio.h>
#include <cmath>
#include <map>
#include <set>
#include <mex.h>
#include "mexFileIO.h"

#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <map>
using namespace std;
#include <mexFileIO.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szParticleF.h>
#include <Feature.h>
#include <DisjointSet.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <Dijkstra.h>
#include <szConvexHull2D.h>

vector<int>
initialClustering(vector<Feature>& P, float space)
{
	vector<Node<int>*> nodes;
	for (int i = 0; i < P.size(); ++i)
	{
		nodes.push_back(makeset(i));
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		for (int j = i + 1; j < nodes.size(); ++j)
		{
			if (Feature::distance(P[i], P[j]) <= space)
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<Node<int>*> c = clusters(nodes);
	vector<int> labels(nodes.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = distance(c.begin(), find(c.begin(), c.end(), findset(nodes[i])));
		labels[i] = k;
	}
	return labels;
}

/*
Divide the feature set (P) into multiple bins, according to the labeling (LABEL).
Small clusters less than (MINSIZE) are not kept.
*/
vector<vector<Feature>> divideClusters(vector<Feature>& P, vector<int>& labels, int minSize)
{
	int numC = 0;
	for (int i = 0; i < labels.size(); ++i)
	{
		numC = Max(labels[i], numC);
	}
	vector<int> counts(numC+1, 0);
	for (int i = 0; i < labels.size(); ++i)
	{
		counts[labels[i]]++;
	}
	vector<vector<Feature>> bins(numC + 1);
	for (int lb = 0; lb <= numC; ++lb)
	{
		if (counts[lb] >= minSize)
		{
			for (int j = 0; j < labels.size(); ++j)
			{
				if (labels[j] == lb)
				{
					bins[lb].push_back(P[j]);
				}
			}
		}
	}
	for (int i = bins.size() - 1; i >= 0; i--)
	{
		if (bins[i].empty())
		{
			bins.erase(bins.begin() + i);
		}
	}
	return bins;
}


/*
Treat the first two elements of a feature as X-Y coordinate, and derive a 2D convex hull from the clustre of features.
*/
vector<Feature> 
convexify(vector<Feature>& cluster)
{
	vector<CParticleF> points;
	for (int i = 0; i < cluster.size(); ++i)
	{
		points.push_back(CParticleF(cluster[i].vals[0], cluster[i].vals[1]));
	}
	vector<CParticleF> h = ConvexHull2D(points);
	vector<Feature> res;

	for (int i = 0; i < h.size(); ++i)
	{
		Feature f(2);
		f.vals[0] = h[i].m_X;
		f.vals[1] = h[i].m_Y;
		res.push_back(f);
	}
	return res;
}
 
/*
From a convex hull, pick a pair that is farthest apart. 
The extremity point can be either of them. 
We pick the first one that appears in the vector.
*/
Feature 
selectExtermity(vector<Feature>& hull)
{
	float dmax = 0;
	Feature extremity;
	for (int i = 0; i < hull.size(); ++i)
	{
		for (int j = i + 1; j < hull.size(); ++j)
		{
			float d = Feature::distance(hull[i], hull[j]);
			if (d > dmax)
			{
				dmax = d;
				extremity = hull[i];
			}
		}
	}
	return extremity;
}



/*
Find shortest paths from each convex vertex to the source.
*/
vector<Vertex<Feature>*>
skeletonTree(vector<Feature>& P, vector<Feature>& hull, Feature& src, float space)
{
	GraphFactory<Feature>& factory = GraphFactory<Feature>::GetInstance();
	vector<Vertex<Feature>*> vertices;
	vector<Vertex<Feature>*> sinks;
	Vertex<Feature>* u = NULL;

	for (int i = 0; i < P.size(); ++i)
	{
		Vertex<Feature>* v = factory.makeVertex(P[i]);
		if (P[i] == src)
		{
			u = v;
		}
		else
		{
			for (int j = 0; j < hull.size(); ++j)
			{
				if (hull[j] == P[i])
				{
					sinks.push_back(v);
					break;
				}
			}
		}
		vertices.push_back(v);
	}
	for (int i = 0; i < vertices.size(); ++i)
	{
		for (int j = i; j < vertices.size(); ++j)
		{
			float d = Feature::distance(P[i], P[j]);
			if (d <= space)
			{
				Edge<Feature>* e1 = factory.makeEdge(vertices[i], vertices[j], d);
				vertices[i]->Add(e1);
				Edge<Feature>* e2 = factory.makeEdge(vertices[j], vertices[i], d);
				vertices[j]->Add(e2);
			}
		}
	}
	Dijkstra(vertices, u);
	set<Vertex<Feature>*> vset;
	for (int i = 0; i < sinks.size(); ++i)
	{
		vector<Vertex<Feature>*> path = tracePath(sinks[i], u);
		for (int j = 0; j < path.size(); ++j)
		{
			vset.insert(path[j]);
		}
	}
	vector<Vertex<Feature>*> tree;
	for (set<Vertex<Feature>*>::iterator it = vset.begin(); it != vset.end(); ++it)
	{
		tree.push_back(*it);
	}
	return tree;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	//printf("%s: This build was compiled at %s %s\n", "OverlapArea", __DATE__, __TIME__);
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: L = EMGClustering(P, K)");
		return;
	}
	//Points
	vector<Feature> P;
	{
		vector<float> P0;
		mxClassID classIdP;
		int ndimP;
		const int* dimsP;
		LoadData(P0, prhs[0], classIdP, ndimP, &dimsP);
		for (int i = 0; i<dimsP[0]; ++i)
		{
			float x = GetData2(P0, i, 0, dimsP[0], dimsP[1], (float)0);
			float y = GetData2(P0, i, 1, dimsP[0], dimsP[1], (float)0);
			Feature f(2);
			f.vals[0] = x;
			f.vals[1] = y;
			P.push_back(f);
		}
	}
	float space = 10.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(space, prhs[1], classMode);
	}
	int minSize = 2;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(minSize, prhs[2], classMode);
	}
	float minDistance = 10.0f;
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(minDistance, prhs[3], classMode);
	}

	vector<int> labels = initialClustering(P, space);
	vector<vector<Feature>> clusters = divideClusters(P, labels, minSize);

	vector<vector<Feature>> clusters2;
	for (int i = 0; i < clusters.size(); ++i)
	{
		vector<Feature> hull = convexify(clusters[i]);
		Feature src = selectExtermity(hull);
		vector<Vertex<Feature>*> t = skeletonTree(clusters[i], hull, src, space);

		vector<Feature> vf;
		for (int k = 0; k < t.size(); ++k)
		{
			vf.push_back(t[k]->key);
		}
		clusters2.push_back(vf);

		/*vector<vector<Vertex<Feature>*>> brs = breakBranhces(t, minSize);
		for (int j = 0; j < brs.size(); ++j)
		{
			vector<Vertex<Feature>*> approx = approximateWithSegments(brs[j], minDistance);
			vector<Feature> vf;
			for (int k = 0; k < approx.size(); ++k)
			{
				vf.push_back(approx[k]->key);
			}
			clusters2.push_back(vf);
		}*/
	}

	if (nlhs >= 1)
	{
		vector<Feature> A;
		vector<Feature> B;
		vector<int> segments;
		for (int i = 0; i < clusters2.size(); ++i)
		{
			for (int j = 0; j < clusters2[i].size()-1; ++j)
			{
				A.push_back(clusters2[i][j]);
				B.push_back(clusters2[i][j+1]);
				segments.push_back(i);
			}
		}
		const int dims[] = { A.size(), 5 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], A[i].vals[0]);
			SetData2(F, i, 1, dims[0], dims[1], A[i].vals[1]);
			SetData2(F, i, 2, dims[0], dims[1], B[i].vals[0]);
			SetData2(F, i, 3, dims[0], dims[1], B[i].vals[1]);
			SetData2(F, i, 4, dims[0], dims[1], (float) segments[i]);
		}
		plhs[0] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		vector<Feature> Q;
		vector<int> segments;
		for (int i = 0; i < clusters2.size(); ++i)
		{
			for (int j = 0; j < clusters2[i].size(); ++j)
			{
				Q.push_back(clusters2[i][j]);
				segments.push_back(i);
			}
		}
		const int dims[] = { Q.size(), 3 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], Q[i].vals[0]);
			SetData2(F, i, 1, dims[0], dims[1], Q[i].vals[1]);
			SetData2(F, i, 2, dims[0], dims[1], (float)segments[i]);
		}
		plhs[1] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		vector<int> labels2;
		for (int i = 0; i < clusters2.size(); ++i)
		{
			for (int j = 0; j < clusters2[i].size(); ++j)
			{
				labels2.push_back(i);
			}
		}
		const int dims[] = { labels2.size(), 1 };
		plhs[2] = StoreData(labels2, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 4) //initial clusters
	{
		vector<Feature> A;
		vector<int> L;
		for (int i = 0; i < clusters.size(); ++i)
		{
			for (int j = 0; j < clusters[i].size() - 1; ++j)
			{
				A.push_back(clusters[i][j]);
				L.push_back(i);
			}
		}
		const int dims[] = { A.size(), 3 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], A[i].vals[0]);
			SetData2(F, i, 1, dims[0], dims[1], A[i].vals[1]);
			SetData2(F, i, 2, dims[0], dims[1], (float)L[i]);
		}
		plhs[3] = StoreData(F, mxSINGLE_CLASS, 2, dims);
	}

	GraphFactory<Feature>::GetInstance().Clean();

	mexUnlock();
}

