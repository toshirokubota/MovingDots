#include <PointLight.h>
#include <szMiscOperations.h>

float
_Distance(PointLight* p, PointLight* q)
{
	float a = 1.0f;
	return sqrt((p->x - q->x)*(p->x - q->x) + (p->y - q->y)*(p->y - q->y) + a * (p->frame - q->frame)*(p->frame - q->frame));
}

/*
Measure how p is compatible to q.
This assumes causality, thus p has to be later frame than q to be compatible.
*/
float
PointLight::_Compatibility(Candidate& a, Candidate& b)
{
	float sgm = 10.0f;
	if (a.tp.q == b.tp.q) return 0.0f; //this happens for points in the first or last frames as well as isolated points.
	float val = 0;
	if (a.tp.p == b.tp.q && a.tp.q == b.tp.r)
	{
		CParticleF a0 = a.tp.evaluate(a.tp.q->frame);
		CParticleF a1 = a.tp.evaluate(a.tp.p->frame);
		CParticleF b0 = b.tp.evaluate(b.tp.q->frame);
		CParticleF b1 = b.tp.evaluate(b.tp.r->frame);
		float d1 = Distance(a0, b1);
		float d2 = Distance(a1, b0);
		val = exp(-(d1*d1 + d2*d2) / (2 * sgm*sgm));
	}
	else if (a.tp.r == b.tp.q && a.tp.q == b.tp.p)
	{
		CParticleF a0 = a.tp.evaluate(a.tp.q->frame);
		CParticleF a1 = a.tp.evaluate(a.tp.r->frame);
		CParticleF b0 = b.tp.evaluate(b.tp.q->frame);
		CParticleF b1 = b.tp.evaluate(b.tp.p->frame);
		float d1 = Distance(a0, b1);
		float d2 = Distance(a1, b0);
		val = exp(-(d1*d1 + d2*d2) / (2 * sgm*sgm));
	}
	else
	{
		val = 0.0f;
	}
	return val;
}

void
PointLight::print(char* tab, char* newl)
{
	if (tab != NULL) printf("%s", tab);
	printf("%d(%3.3f,%3.3f,%f) [", id, x, y, frame);
	if (newl != NULL) printf("%s", newl);
}

void
PointLight::initializeProb()
{
	int n = candidates.size();
	for (int i = 0; i < candidates.size(); ++i)
	{
		candidates[i].prob = 1.0 / (n + 1.0);
	}
}

PointLight::PointLightTriple 
PointLight::winner()
{ 
	float maxP = 0;
	float totalP = 0;
	int idx = -1;
	for (int i = 0; i < candidates.size(); ++i)
	{
		if (candidates[i].prob > maxP)
		{
			maxP = candidates[i].prob;
			idx = i;
		}
		totalP += candidates[i].prob;
	}
	float eta = 1.0 - totalP;
	if (eta > maxP)
	{
		return PointLightTriple(NULL, this, NULL);
	}
	else
	{
		return candidates[idx].tp;
	}
}

void
PointLight::updateFitness()
{
	for (int i = 0; i < candidates.size(); ++i)
	{
		float s = 0;
		for (int j = 0; j < candidates[i].tp.p->candidates.size(); ++j)
		{
			float val = _Compatibility(candidates[i], candidates[i].tp.p->candidates[j]) * candidates[i].tp.p->candidates[j].prob;
			s += val;
			/*if (val > 0)
			{
				if ((id == 39 && (i == 6 || i == 8)))
				{
					printf("(%d %d %d)- (%d %d %d) %f %f\n",
						candidates[i].tp.p->id, candidates[i].tp.q->id, candidates[i].tp.r->id,
						candidates[i].tp.p->candidates[j].tp.p->id, candidates[i].tp.p->candidates[j].tp.q->id, candidates[i].tp.p->candidates[j].tp.r->id,
						val, s);
				}
			}*/
			//s = Max(s, _Compatibility(candidates[i], candidates[i].tp.p->candidates[j]) * candidates[i].tp.p->candidates[j].prob);
		}
		for (int j = 0; j < candidates[i].tp.r->candidates.size(); ++j)
		{
			float val = _Compatibility(candidates[i], candidates[i].tp.r->candidates[j]) * candidates[i].tp.r->candidates[j].prob;
			s += val;
			/*if (val > 0)
			{
				if ((id == 39 && (i == 6 || i == 8)))
				{
					printf("(%d %d %d)- (%d %d %d) %f %f\n", 
						candidates[i].tp.p->id, candidates[i].tp.q->id, candidates[i].tp.r->id,
						candidates[i].tp.r->candidates[j].tp.p->id, candidates[i].tp.r->candidates[j].tp.q->id, candidates[i].tp.r->candidates[j].tp.r->id,
						val, s);
				}
			}*/
			//s = Max(s, _Compatibility(candidates[i], candidates[i].tp.r->candidates[j]) * candidates[i].tp.r->candidates[j].prob);
		}
		candidates[i].comp = s;
		/*if (id == 99 || id == 100 || id == 40 || id == 39 || id == 41 || id == 101)
		{
			if (candidates[i].comp > 0.1)
			{
				printf("F %d) %d (%d %d %d) => %f (%f)\n", id, i, candidates[i].tp.p->id, candidates[i].tp.q->id, candidates[i].tp.r->id,
					candidates[i].comp, candidates[i].prob);
			}
		}*/
	}
}

void
PointLight::updateProb()
{
	float totalP = 0;
	float total = 0;
	for (int i = 0; i < candidates.size(); ++i)
	{
		total += candidates[i].comp * candidates[i].prob;
		totalP += candidates[i].prob;
	}
	float eta = 0; // totalP >= 1.0 ? 0.0f : (1.0 - totalP) * Threshold;
	for (int i = 0; i < candidates.size(); ++i)
	{
		candidates[i].prob = candidates[i].comp * candidates[i].prob / (total + eta);
		/*if (id == 99 || id == 100 || id == 40 || id == 39 || id==41 || id==101)
		{
			if (candidates[i].prob > 0.1)
			{
				printf("P %d) %d (%d %d %d) => %f (%f)\n", id, i, candidates[i].tp.p->id, candidates[i].tp.q->id, candidates[i].tp.r->id, 
					candidates[i].comp, candidates[i].prob);
			}
		}*/
	}
}
