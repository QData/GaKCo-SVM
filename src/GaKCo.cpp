// GaKCo : Fast Gapped k-mer string Kernel using Counting
// Code Contibution by:
//Ritambhara Singh <rs3zz@virginia.edu>
//Kamran Kowsari <kk7nc@virginia.edu >
//Arshdeep Sekhon <as5cu@virginia.edu >


// This file contains Main Code


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "shared.h"
#include "shared.cpp"
#include <assert.h>
#include <omp.h>
#include <thread>
#include <iostream>
#include <random>
#include <ctime>
#include "readInput.cpp"
#include <future>
#define ARG_REQUIRED 7


//extract g-mers from input sequences
Features *extractFeatures(int **S, int *len, int nStr, int g)
{
	 int i, j, j1;
	 int n, sumLen, nfeat, addr;
	int *group;
	int *features;

	int *s;
	 int c;
	Features *F;

	nfeat = 0;
	sumLen = 0;

	for (i = 0; i < nStr; ++i)
	{
		sumLen += len[i];
		nfeat += (len[i] >= g) ? (len[i] - g + 1) : 0;
	}

	printf("numF=%d, sumLen=%d\n", nfeat, sumLen); 
	group = (int *)malloc(nfeat * sizeof(int));
	features = (int *)malloc(nfeat*g * sizeof(int *));

	c = 0;

	for (i = 0; i < nStr; ++i)
	{
		s = S[i];



		for (j = 0; j < len[i] - g + 1; ++j)
		{

			for (j1 = 0; j1 <g; ++j1)
			{
				features[c + j1*nfeat] = s[j + j1];
				
			}
			group[c] = i;
			c++;
		}
	}
	if (nfeat != c)
		printf("Something is wrong...\n");

	F = (Features *)malloc(sizeof(Features));
	(*F).features = features;
	(*F).group = group;
	(*F).n = nfeat;
	return F;
}


int help()
{
	printf("Usage: trainKernel <Sequence-file> <Dictionary> <Label_output> <g> <k>  <Kernel-file> <setParallel>\n");
	printf("\t Sequence-file : file with sequence data\n");
	printf("\t g : length of gapped instance, >0 and < 20 \n");
	printf("\t k : length of k-mer, < g");
	printf("\t g-k < 20");
	printf("setParallel : 1 for multithread");
	printf("\t IMPORTANT: sequence elements must be\n\tin the range [0,AlphabetSize - 1].\n");
	return 1;
}

int errorID1()
{
	printf("Error: g >= Shortest sequence in the input file!\n");
	return 1;
}



// Build cumulative mismatch profile for each M

void main_loop_kernel(int * elems,Features * features ,unsigned int *Ksfinal,int * cnt_k,   int *feat,int g, int k, int na,int nfeat,int nStr,int i)
{

	unsigned long int c = 0;
	int num_comb;
	Combinations * combinations = (Combinations *)malloc(sizeof(Combinations));
        unsigned int *Ks = (unsigned int *)malloc(nStr*nStr * sizeof(unsigned int));
	unsigned int *sortIdx = (unsigned int *)malloc(nfeat * sizeof(unsigned int));
	unsigned int *features_srt  = (unsigned int *)malloc(nfeat*g * sizeof(unsigned int *));
	unsigned int *group_srt = (unsigned int *)malloc(nfeat * sizeof(unsigned int));
	unsigned int *cnt_comb = (unsigned int *)malloc(2 * sizeof(unsigned int));
	unsigned int *feat1 = (unsigned int *)malloc(nfeat*g * sizeof(unsigned int));
	
	int *pos = (int *)malloc(nfeat * sizeof(int));
	memset(pos, 0, sizeof(int) * nfeat);
	c =i*(nStr*nStr);


			(*combinations).n = g;
			(*combinations).k = g - i;

			(*combinations).num_comb = nchoosek(g, g - i);

			// number of possible positions
			num_comb = nchoosek(g, g - i);
 
			unsigned int  *out = (unsigned int *)malloc((g - i)*num_comb * sizeof(unsigned int));
			unsigned int  *cnt_m = (unsigned int *)malloc(g * sizeof(unsigned int));

			cnt_comb[0] = 0;

			getCombinations(elems,(*combinations).n, (*combinations).k, pos, 0, 0, cnt_comb, out, num_comb);
			cnt_m[i] = cnt_comb[0];

			cnt_comb[0] += ((*combinations).k*num_comb);
				for ( int j = 0; j < num_comb; ++j)
				{
					//remove i positions
						for ( int j1 = 0; j1 < nfeat; ++j1)
						{
								for ( int j2 = 0; j2 < g - i; ++j2)
								{
									feat1[j1 + j2*nfeat] = feat[j1 + (out[(cnt_m[i] - num_comb + j) + j2*num_comb])*nfeat];
									
								}
							
							
						}
						//sort the g-mers
					       cntsrtna(sortIdx,feat1, g - i, nfeat, na);
					    
					       
						for ( int j1 = 0; j1 < nfeat; ++j1)
						{

								for ( int j2 = 0; j2 < g - i; ++j2)
								{
									features_srt[j1 + j2*nfeat] = feat1[(sortIdx[j1]) + j2*nfeat];
									
								}
							group_srt[j1] = (*features).group[sortIdx[j1]];
						}
					//update cumulative mismatch profile
					countAndUpdate(Ks, features_srt, group_srt, g - i, nfeat, nStr);
					
						for ( int j1 = 0; j1 < nStr; ++j1)
						{
							for ( int j2 = j1; j2 < nStr; ++j2)
								
							{
							      if(j1!=j2)
								Ksfinal[(c + j1) + j2*nStr] +=  Ks[j1 + j2*nStr];
								Ksfinal[c +(j1)*nStr + j2] +=  Ks[j1 + j2*nStr];
								
							}
						}
					
					

					
				}
				
			free(cnt_m);
			free(out);
			cnt_k[i] = c;
			printf("iter:%d\n", i); 
			free(Ks);
			free(sortIdx);
			free(features_srt);
			free(group_srt);
			free(feat1);
			free(cnt_comb);
			free(pos);
}



//Main function 

int main(int argc, char *argv[])
{
	 	
	char filename[100],filename_label[100],Dicfilename[100], opfilename[100];
	int *label;
	int k, num_max_mismatches;
	int m, g;
	int na;
	unsigned int addr;
	long int nStr, num_comb, value;
	
	int nfeat;
	double *K;
	unsigned int *nchoosekmat, *Ks, *Ksfinal, *Ksfinalmat;
	int *len;
	int **S;
	unsigned int *sortIdx;
	int *feat;
	unsigned int  *out, *out_temp, *resgroup;
	int *elems, *cnt_k;
	unsigned int *cnt_m;
	long int maxIdx, maxlen, minlen;
	char isVerbose;
	Features *features;
	int parallel_;
	isVerbose = 0;
	if (argc != ARG_REQUIRED + 1)
	{
		return help();
	}
	strcpy(filename ,argv[1]);
	strcpy(Dicfilename ,argv[2]);
	strcpy(filename_label ,argv[3]);
	g = atoi(argv[4]);
	k = atoi(argv[5]);
	strcpy(opfilename, argv[6]);
	parallel_=atoi(argv[7]);

	
	
	label = (int *)malloc(STRMAXLEN*sizeof(int));
	len = (int *)malloc(STRMAXLEN * sizeof(int));
	assert(len != 0);  
	maxlen = 0;
	minlen = MAXNSTR;
	nStr = MAXNSTR;
	
	// Read input file
	
	printf("Input file : %s\n", filename);
	S = Readinput_(filename,Dicfilename,label,len, &nStr, &maxlen, &minlen,&na);
	
	
	
	if (k <= 0 || g <= k || g>20 || g - k>20 || na <= 0)
		return help();
	if (maxlen != minlen)
		printf("Read %ld strings of max length = %ld and min length=%ld\n", nStr, maxlen, minlen);
	else
		printf("Read %ld strings of length = %ld\n", nStr, maxlen);

	if (g > minlen)
		return errorID1();

	/* Precompute weights hm.*/

	int w[g - k];
	printf("Weights (hm):");
		for (int i = 0; i <= g - k; i++)
		{
			w[i] = nchoosek(g - i, k);
			printf("%d ", w[i]);
		}
	
	printf("\n");

	/*Extract g-mers.*/
	features = extractFeatures(S, len, nStr, g);
	
	nfeat = (*features).n;
	feat = (*features).features;
	printf("(%d,%d): %d features\n", g, k, nfeat); 


	/*Compute gapped kernel.*/
	K = (double *)malloc(nStr*nStr * sizeof(double));

	

	addr = ((g - k) + 1)*nStr*nStr;
	
	Ksfinal = (unsigned int *)malloc(addr * sizeof(unsigned int));
	
	memset(Ksfinal, 0, sizeof(unsigned int) * addr);
	

	
	elems = (int *)malloc(g * sizeof(int));

	
	cnt_k = (int *)malloc(nfeat * sizeof(int));
		for (int i = 0; i < g; ++i)
		{
			elems[i] = i;
		}

std::vector<std::thread> th;
for ( int i = 0 ; i <= g-k; i++)
	{
	  if(parallel_)
	  {
	    th.push_back(std::thread(&main_loop_kernel,elems,features ,Ksfinal,cnt_k,feat,g, k, na,nfeat,nStr,i));
	  }else
	  {
	  main_loop_kernel(elems,features ,Ksfinal,cnt_k,feat,g, k, na,nfeat,nStr,i);
	  }
	}
	
	nchoosekmat = (unsigned int *)malloc(g*g * sizeof(unsigned int));
	memset(nchoosekmat, 0, sizeof(unsigned int) * g * g);
	

	
	
	

	
	for ( int i = g; i >= 0; --i)
		{
				for ( int j = 1; j <= i; ++j)
				{
				  
					nchoosekmat[(i - 1) + (j - 1)*g] = nchoosek(i, j);
				      
				}
			
		}





	 int c1 = 0,
	c2 = 0;

	num_max_mismatches = g - k;

	//join parallel threads
if(parallel_)
{
  for(auto &t : th)
  {
	  t.join();
  }
}
//get exact mismatch profile
		for ( int i = 1; i <= num_max_mismatches; ++i)
		{
			c1 = cnt_k[i];

			for ( int j = 0; j <= i - 1; ++j)
			{
				c2 = cnt_k[j];
				for ( int j1 = 0; j1 < nStr; ++j1)
				{
					value = 0;
					int x = 0;

					for ( int j2 = 0; j2 < nStr; ++j2)
					{
						Ksfinal[(c1 + j1) + j2*nStr] -=  nchoosekmat[(g - j - 1) + (i - j - 1)*g] * Ksfinal[(c2 + j1) + j2*nStr];
						
					}
				}
			}
		}
	
		for ( int i = 0; i <= g - k; i++)
		{
			c1 = cnt_k[i];
			for ( int j1 = 0; j1 < nStr; ++j1)
			{
				for ( int j2 = 0; j2 < nStr; ++j2)
				{
				 K[j1 + j2*nStr] += w[i] * Ksfinal[(c1 + j1) + j2*nStr];

				}
				
			}
		}
	/*Normalize kernel values and write into a file*/
	FILE *kernelfile;
	FILE *labelfile;
	kernelfile = fopen(opfilename, "w");
	labelfile = fopen(filename_label, "w");


		for (int i = 0; i < nStr; ++i)
		{	
			for (int j = 0; j < nStr; ++j)
			{
				fprintf(kernelfile, "%d:%e ", j + 1, K[i + j*nStr] / sqrt(K[i + i*nStr] * K[j + j*nStr]));
			}
			fprintf(labelfile, "%d ", label[i]);
			fprintf(labelfile, "\n");
			fprintf(kernelfile, "\n");
		}
		
	fclose(kernelfile);
	fclose(labelfile);
	free(cnt_k);
	free(Ksfinal);
	free(label);
	free(K);
	free(nchoosekmat);
	free(feat);
	free(elems);
	return 0;
}
