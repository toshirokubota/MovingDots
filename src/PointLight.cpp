#include <PointLight.h>

float
_Distance(PointLight* p, PointLight* q)
{
	float a = 1.0f;
	return sqrt((p->x - q->x)*(p->x - q->x) + (p->y - q->y)*(p->y - q->y) + a * (p->frame - q->frame)*(p->frame - q->frame));
}

/*
Compute the average location of p at the frame of q.
Then find the distance between the two locations at the frame.
*/
float
_PredictedDistance(PointLight* p, PointLight* q)
{
	float vx = 0, vy = 0;
	for (int i = 0; i < p->prevC.size(); ++i)
	{
		vx += p->prevP[i] * (p->x - p->prevC[i]->x) / (p->frame - p->prevC[i]->frame);
		vy += p->prevP[i] * (p->y - p->prevC[i]->y) / (p->frame - p->prevC[i]->frame);
	}
	float xp = p->x + vx * (q->frame - p->frame);
	float yp = p->y + vy * (q->frame - p->frame);
	return sqrt((xp - q->x)*(xp - q->x) + (yp - q->y)*(yp - q->y));
}

float
_Continuity(PointLight* p, PointLight* q, int m, int n, float sgm)
{
	vector<float> pvals;
	vector<PointLight*> cand;
	if (p->frame < q->frame)
	{
		pvals = p->prevP;
		cand = p->prevC;
	}
	else if(p->frame > q->frame)
	{
		pvals = p->nextP;
		cand = p->nextC;
	}
	else {
		return 0;
	}

	//find the average landing point given the candidates
	float nx = 0, ny = 0, nf = q->frame - p->frame;
	for (int i = 0; i < pvals.size(); ++i)
	{
		PointLight* r = i < cand.size() ? cand[i] : p;
		float df = r->frame - p->frame;
		if (Abs(df) > 0)
		{
			float vx = (r->x - p->x) / df;
			float vy = (r->y - p->y) / df;
			nx += pvals[i] * vx * nf;
			ny += pvals[i] * vy * nf;
		}
	}
	float px = p->x + nx;  //predicted x
	float py = p->y + ny;  //predicted y
	float ee = (q->x - px) * (q->x - px) + (q->y - py)*(q->y - py);
	return exp(-ee / (2.*sgm*sgm));
}

float
_Compatibility(PointLight* p, PointLight* q, int m, int n, float sgm)
{
	float c = (_Continuity(p, q, m, n, sgm) + _Continuity(q, p, n, m, sgm)) / 2.0f;
	float d = _Distance(p, q);
	float w = exp(-d*d / (2 * sgm*sgm));
	float ee = w * d + (1.0 - w)* c;
	return exp(-ee * ee/ (2.*sgm*sgm));
}

PointLight*
PointLight::nextWinner()
{
	PointLight* winner = NULL;
	if (nextC.size() > 0)
	{
		float maxp = nextP[nextC.size()];
		for (int i = 0; i < nextC.size(); ++i)
		{
			if (nextP[i] > maxp)
			{
				winner = nextC[i];
				maxp = nextP[i];
			}
		}
	}
	return winner;
}

PointLight*
PointLight::prevWinner()
{
	PointLight* winner = NULL;
	if (prevC.size() > 0)
	{
		float maxp = prevP[prevC.size()];
		for (int i = 0; i < prevC.size(); ++i)
		{
			if (prevP[i] > maxp)
			{
				winner = prevC[i];
				maxp = prevP[i];
			}
		}
	}
	return winner;
}


void
PointLight::print(char* tab, char* newl)
{
	if (tab != NULL) printf("%s", tab);
	printf("%d(%3.3f,%3.3f,%f) [", id, x, y, frame);
	for (int i = 0; i < prevC.size(); ++i)
	{
		printf("%3.3f(%d) ", prevP[i], prevC[i]->id);
	}
	printf("%3.3f], [", prevP[prevC.size()]);
	for (int i = 0; i < nextC.size(); ++i)
	{
		printf("%3.3f(%d) ", nextP[i], nextC[i]->id);
	}
	printf("%3.3f]", nextP[nextC.size()]);
	if (newl != NULL) printf("%s", newl);
}

/*
Using previous and future frames, find likely candidates for this point in both previous and future frames.
This is used to initialize points at the beginning of a video where little velocity informaiton is available.
*/
void
PointLight::velocityFreeInitialization(vector<PointLight*>& points)
{
	vector<pair<float, PointLight*>> prevPairs;
	vector<pair<float, PointLight*>> nextPairs;
	int  range = 10;
	for (int i = 0; i < points.size(); ++i)
	{
		if (points[i] == this) continue;
		float df = this->frame - points[i]->frame;
		if (Abs(df) <= range)
		{
			if (df > 0)
			{
				prevPairs.push_back(pair<float, PointLight*>(_Distance(this, points[i]), points[i]));
			}
			else if (df < 0)
			{
				nextPairs.push_back(pair<float, PointLight*>(_Distance(this, points[i]), points[i]));
			}
		}
	}
	sort(prevPairs.begin(), prevPairs.end());
	sort(nextPairs.begin(), nextPairs.end());
	for (int i = 0; i < Min(prevPairs.size(), NumCandidates); ++i)
	{
		PointLight* p = prevPairs[i].second;
		prevC.push_back(p);
	}
	for (int i = 0; i < Min(nextPairs.size(), NumCandidates); ++i)
	{
		PointLight* p = nextPairs[i].second;
		nextC.push_back(p);
	}

	prevP = vector<float>(prevC.size() + 1, 1.0 / (prevC.size() + 1));
	nextP = vector<float>(nextC.size() + 1, 1.0 / (nextC.size() + 1));
	prevP0 = prevP;
	nextP0 = nextP;
}

/*
Using previous frames, find likely candidates for this point in the previous frames.
*/
void
PointLight::velocityDrivenInitialization(vector<PointLight*>& points)
{
	vector<pair<float, PointLight*>> prevPairs;
	int  range = 10;
	for (int i = 0; i < points.size(); ++i)
	{
		if (points[i] == this) continue;
		float df = this->frame - points[i]->frame;
		if (Abs(df) <= range)
		{
			if (df > 0)
			{
				prevPairs.push_back(pair<float, PointLight*>(_PredictedDistance(this, points[i]), points[i]));
			}
		}
	}
	sort(prevPairs.begin(), prevPairs.end());
	for (int i = 0; i < Min(prevPairs.size(), NumCandidates); ++i)
	{
		PointLight* p = prevPairs[i].second;
		prevC.push_back(p);
		p->nextC.push_back(this);
	}
	prevP = vector<float>(prevC.size() + 1, 1.0 / (prevC.size() + 1));
	prevP0 = prevP;

	for (int i = 0; i < prevC.size(); ++i)
	{
		prevC[i]->nextC.push_back(this);
		prevC[i]->nextP = vector<float>(prevC[i]->nextC.size() + 1, 1.0 / (prevC[i]->nextC.size() + 1));
		prevC[i]->nextP0 = prevC[i]->nextP;
	}
}

void
PointLight::update0(float sigma, float thres, float rate)
{
	{
		int numc1 = this->prevC.size();
		if (numc1)
		{
			float total_sum = 0;
			vector<float> comp(numc1, 0.0f);
			for (int m = 0; m < numc1; ++m)
			{
				PointLight* p = this->prevC[m];
				for (int m2 = 0; m2 < p->nextC.size(); ++m2)
				{
					if (p->nextC[m2] == this)
					{
						comp[m] += p->nextP[m2] * _Compatibility(this, p, m, m2, sigma);
						break;
					}
				}
				total_sum += prevP[m] * comp[m];
			}
			if (total_sum > 1.0e-5)
			{
				for (int m = 0; m < numc1; ++m)
				{
					prevP0[m] = prevP[m] * comp[m] / (total_sum + prevP[numc1] * thres);
				}
				prevP0[numc1] = prevP[numc1] * thres / (total_sum + prevP[numc1] * thres);
			}
			else
			{
				for (int m = 0; m < numc1; ++m)
				{
					this->prevP0[m] = 0.0f;
				}
				this->prevP0[numc1] = 1.0f;
			}
		}
	}
	
	{
		int numc1 = this->nextC.size();
		if (numc1 > 0)
		{
			float total_sum = 0;
			vector<float> comp(numc1, 0.0f);
			for (int m = 0; m < numc1; ++m)
			{
				PointLight* p = this->nextC[m];
				for (int m2 = 0; m2 < p->prevC.size(); ++m2)
				{
					if (p->prevC[m2] == this)
					{
						comp[m] += p->prevP[m2] * _Compatibility(this, p, m, m2, sigma);
						break;
					}
				}
				total_sum += nextP[m] * comp[m];
			}
			if (total_sum > 1.0e-5)
			{
				for (int m = 0; m < numc1; ++m)
				{
					nextP0[m] = nextP[m] * comp[m] / (total_sum + nextP[numc1] * thres);
				}
				nextP0[numc1] = nextP[numc1] * thres / (total_sum + nextP[numc1] * thres);
			}
			else
			{
				for (int m = 0; m < numc1; ++m)
				{
					this->nextP0[m] = 0.0f;
				}
				this->nextP0[numc1] = 1.0f;
			}
		}
	}
}


