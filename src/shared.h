#pragma once



typedef struct Feature
{
	int *features;
	int *group;
	int n;
} Features;

typedef struct Combinations
{
	int n;
	int k;
	double num_comb;
	int *comb;
} Combinations;

typedef struct Dict
{
	int index;
	char *word;
} Dict;





char *trimwhitespace(char *s);
void cntsrtna(unsigned int *out,unsigned int *sx, int k, int r, int na);
void countAndUpdate(unsigned int *outK, unsigned int *sx, unsigned int *g, int k, int r, int nStr);
double nchoosek(double n, double k);
void getCombinations(const int *elems, unsigned int n, unsigned int k, int *pos, unsigned int depth, unsigned int margin, unsigned int *cnt_comb, unsigned int *out, int num_comb);

