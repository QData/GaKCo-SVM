// GaKCo : Fast Gapped k-mer string Kernel using Counting
// Code Contibution by:
//Ritambhara Singh <rs3zz@virginia.edu>
//Kamran Kowsari <kk7nc@virginia.edu>
//Arshdeep Sekhon <as5cu@virginia.edu>
//Derrick Blakely <dcb7xz@virginia.edu>

// This file contains Main Code


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "shared.h"
#include "shared.cpp"
#include <assert.h>
#include <thread>
#include <iostream>
#include <random>
#include <ctime>
#include "readInput.cpp"
#include <future>
#include <unistd.h>

//extract g-mers from input sequences
Features *extractFeatures(int **S, int *len, int nStr, int g) {
	int i, j, j1;
	int n, nfeat, addr;
	int *group;
	int *features;
	int *s;
	int c;
	Features *F;

	nfeat = 0;
	for (i = 0; i < nStr; ++i) {
		nfeat += (len[i] >= g) ? (len[i] - g + 1) : 0;
	}

	printf("Number of gmers: %d\n", nfeat); 
	group = (int *)malloc(nfeat * sizeof(int));
	features = (int *)malloc(nfeat*g * sizeof(int *));
	c = 0;

	for (i = 0; i < nStr; ++i) {
		s = S[i];
		for (j = 0; j < len[i] - g + 1; ++j) {
			for (j1 = 0; j1 < g; ++j1) {
				features[c + j1*nfeat] = s[j + j1];
			}
			group[c] = i;
			c++;
		}
	}
	if (nfeat != c) {
		printf("Something is wrong...\n");
	}

	F = (Features *)malloc(sizeof(Features));
	(*F).features = features;
	(*F).group = group;
	(*F).n = nfeat;
	return F;
}


int help() {
	printf("\nUsage: gakco [options] <sequenceFile> <dictionaryFile> <labelsFile> <kernelFile>\n");
	printf("\t g : length of gapped instance. Constraints: 0 < g < 20\n");
	printf("\t k : length of k-mer. Constraints: k < g\n");
	printf("\t n : (optional) maximum number of examples in the data set. Default: 15000\n");
	printf("\t p : (optional) parallel. Set to 1 to using multithreading; else set to 0. Default: 1\n");
	printf("\t sequenceFile : set of sequences. Format must be as follows:\n");
	printf("\t\t>label\n\t\tsequence\n\t\tetc.\n");
	printf("\t\tWhere the label is a class label (-1 or 0 for negative sequences, 1 for positive sequences).\n");
	printf("\t\tSequences must be contained in a single line. See \"Sequence File Format\" section of readme for more details.\n");
	printf("\t dictionaryFile : file containing the alphabet of characters that appear in the sequences (simple text file)\n");
	printf("\t labelsFile : file to place labels from the examples (simple text file)\n");
	printf("\t kernelFile : name of the file to write the kernel that will be computed by GaKCo\n");
	printf("\n");
	printf("\t IMPORTANT: \n");
	printf("\t\t sequence elements must be in the range [0,AlphabetSize - 1].\n");
	printf("\t\t g - k should be less than 20\n");
	printf("\nExample usage: ./GaKCo -g 7 -k 5 -n 15000 -p 1 myProteinFile.fasta myProteinDictionary.txt labelOutputFile.txt kernelOutputFile.txt\n\n");

	return 1;
}

int errorID1() {
	printf("Error: g >= Shortest sequence in the input file!\n");
	return 1;
}

// Build cumulative mismatch profile for each M

void main_loop_kernel(int * elems,Features * features ,unsigned int *Ksfinal,int * cnt_k,   int *feat,int g, int dictionarySize, int nfeat,int nStr,int i) {
	unsigned long int c = 0;
	int num_comb;
	Combinations * combinations = (Combinations *)malloc(sizeof(Combinations));
	unsigned int *Ks = (unsigned int *)malloc(nStr*nStr * sizeof(unsigned int));
	unsigned int *sortIdx = (unsigned int *)malloc(nfeat * sizeof(unsigned int));
	unsigned int *features_srt  = (unsigned int *)malloc(nfeat*g * sizeof(unsigned int *));
	unsigned int *group_srt = (unsigned int *)malloc(nfeat * sizeof(unsigned int));
	unsigned int *cnt_comb = (unsigned int *)malloc(2 * sizeof(unsigned int));
	unsigned int *feat1 = (unsigned int *)malloc(nfeat*g * sizeof(unsigned int));
	
	int *pos = (int *)malloc(g * sizeof(int));
	memset(pos, 0, sizeof(int) * g);

	c = i*(nStr*nStr);
	(*combinations).n = g;
	(*combinations).k = g - i;
	int k = g - i;

	(*combinations).num_comb = nchoosek(g, k);

	// number of possible positions
	num_comb = nchoosek(g, k);
	unsigned int  *out = (unsigned int *)malloc(k*num_comb * sizeof(unsigned int));
	unsigned int  *cnt_m = (unsigned int *)malloc(g * sizeof(unsigned int));
	cnt_comb[0] = 0;
	getCombinations(elems,(*combinations).n, (*combinations).k, pos, 0, 0, cnt_comb, out, num_comb);
	cnt_m[i] = cnt_comb[0];

	cnt_comb[0] += ((*combinations).k*num_comb);
	for ( int j = 0; j < num_comb; ++j) {
		//remove i positions 
		for ( int j1 = 0; j1 < nfeat; ++j1) {
			for ( int j2 = 0; j2 < k; ++j2) {
				unsigned int out_val = out[(cnt_m[i] - num_comb + j) + j2*num_comb];
				feat1[j1 + j2*nfeat] = feat[j1 + out_val*nfeat];
			}
		}
		//sort the g-mers
		cntsrtna(sortIdx,feat1, k, nfeat, dictionarySize);    

		for ( int j1 = 0; j1 < nfeat; ++j1) {
			for ( int j2 = 0; j2 < k; ++j2) {
				features_srt[j1 + j2*nfeat] = feat1[(sortIdx[j1]) + j2*nfeat];
			}
			group_srt[j1] = (*features).group[sortIdx[j1]];
		}
		//update cumulative mismatch profile
		countAndUpdate(Ks, features_srt, group_srt, k, nfeat, nStr);
		
		for ( int j1 = 0; j1 < nStr; ++j1) {
			for ( int j2 = j1; j2 < nStr; ++j2) {
				if (j1 != j2) {
					Ksfinal[(c + j1) + j2*nStr] += Ks[j1 + j2*nStr];
				}
				Ksfinal[c + j1*nStr + j2] += Ks[j1 + j2*nStr];
			}
		}
	}
	
	free(cnt_m);
	free(out);
	cnt_k[i] = c;
	printf("Finished handling mismatch level of m = %d\n", i); 
	free(Ks);
	free(sortIdx);
	free(features_srt);
	free(group_srt);
	free(feat1);
	free(cnt_comb);
	free(pos);
}

//Main function 
int main(int argc, char *argv[]) {
	//Get the g, k, nStr, and parallel_ values from the command line 
	int g = -1;
	int k = -1;
	int parallel_ = 1; //use multithreading by default
	long int maxNumStr = 15000; //Use max string length of 15000 by default

	int c;
	while ((c = getopt(argc, argv, "g:k:n:p:")) != -1) {
		switch (c) {
			case 'g':
				g = atoi(optarg);
				break;
			case 'k':
				k = atoi(optarg);
				break;
			case 'n':
				maxNumStr = atoi(optarg);
				break;
			case 'p':
				parallel_ = atoi(optarg);
				if (parallel_ != 0 && parallel_ != 1) {
					return help();
				}
				break;
		}
	}
	if (g == -1) {
		printf("Must provide a value for the g parameter\n");
		return help();
	}
	if (k == -1) {
		printf("Must provide a value for the k parameter\n");
		return help();
	}

	int argNum = optind;
	if (argc - argNum != 4) {
		printf("Incorrect number of arguments were provided.\n");
		return help();
	}

	//Get the names of the sequenceFile, dictionaryFile, labelsFile, and kernelFile from the command line
	char filename[100], filename_label[100], dictFileName[100], opfilename[100];
	strcpy(filename, argv[argNum++]);
	strcpy(dictFileName, argv[argNum++]);
	strcpy(filename_label, argv[argNum++]);
	strcpy(opfilename, argv[argNum]);

	int *seqLabels;
	int num_max_mismatches, dictionarySize;
	unsigned int addr;
	long int num_comb;
	int nfeat;
	double *K;
	unsigned int *nchoosekmat, *Ks, *Ksfinal, *Ksfinalmat;
	int *seqLengths;
	int **S;
	unsigned int *sortIdx;
	int *feat;
	unsigned int  *out, *out_temp, *resgroup;
	int *elems, *cnt_k;
	unsigned int *cnt_m;
	long int maxIdx, maxlen, minlen;
	char isVerbose;
	Features *features;
	
	isVerbose = 0;
	
	seqLabels = (int *) malloc(maxNumStr * sizeof(int));
	seqLengths = (int *) malloc(maxNumStr * sizeof(int));
	assert(seqLengths != 0);  
	maxlen = 0;
	minlen = STRMAXLEN;
	long int nStr = maxNumStr;
	
	//Read the sequence file
	
	printf("Input file : %s\n", filename);
	S = Readinput_(filename, dictFileName, seqLabels, seqLengths, &nStr, &maxlen, &minlen, &dictionarySize, maxNumStr);
	
	if (k <= 0 || g <= k || g > 20 || g - k > 20 || dictionarySize <= 0) {
		return help();
	}
	if (maxlen != minlen) {
		printf("Read %ld strings of max length = %ld and min length=%ld\n", nStr, maxlen, minlen);
	} else {
		printf("Read %ld strings of length = %ld\n", nStr, maxlen);
	}
	if (g > minlen) {
		return errorID1();
	}

	num_max_mismatches = g - k;

	/* Precompute weights hm.*/
	int w[g - k + 1];
	for (int i = 0; i <= num_max_mismatches; i++) {
		w[i] = nchoosek(g - i, k);
	}

	/*Extract g-mers.*/
	features = extractFeatures(S, seqLengths, nStr, g);
	
	nfeat = (*features).n;
	feat = (*features).features;

	/*Compute gapped kernel.*/
	K = (double *)malloc(nStr*nStr * sizeof(double));
	memset(K, 0, sizeof(double *) * nStr * nStr);

	addr = ((g - k) + 1) * nStr * nStr;

	Ksfinal = (unsigned int *)malloc(addr * sizeof(unsigned int));
	
	memset(Ksfinal, 0, sizeof(unsigned int) * addr);
		
	elems = (int *) malloc(g * sizeof(int));

	cnt_k = (int *) malloc((g - k + 1) * sizeof(int));
	for (int i = 0; i < g; ++i) {
		elems[i] = i;
	}

	std::vector<std::thread> th;
	for ( int i = 0 ; i <= g-k; i++) {
		if(parallel_) {
			th.push_back(std::thread(&main_loop_kernel,elems,features ,Ksfinal,cnt_k,feat,g, dictionarySize, nfeat,nStr,i));
		} else {
			main_loop_kernel(elems,features, Ksfinal, cnt_k, feat, g, dictionarySize, nfeat, nStr, i);
		}
	}
	
	nchoosekmat = (unsigned int *)malloc(g*g * sizeof(unsigned int));
	memset(nchoosekmat, 0, sizeof(unsigned int) * g * g);
	
	for ( int i = g; i >= 0; --i) {
		for ( int j = 1; j <= i; ++j) {
			nchoosekmat[(i - 1) + (j - 1)*g] = nchoosek(i, j);
		}
	}

	int c1 = 0;
	int c2 = 0;

	//join parallel threads
	if(parallel_) {
		for(auto &t : th) {
			t.join();
		}
	}
	//get exact mismatch profile
	for ( int i = 1; i <= num_max_mismatches; ++i) {
		c1 = cnt_k[i];
		for ( int j = 0; j <= i - 1; ++j) {
			c2 = cnt_k[j];
			for ( int j1 = 0; j1 < nStr; ++j1) {
				for ( int j2 = 0; j2 < nStr; ++j2) {
					Ksfinal[(c1 + j1) + j2*nStr] -=  nchoosekmat[(g - j - 1) + (i - j - 1)*g] * Ksfinal[(c2 + j1) + j2*nStr];
				}
			}
		}
	}
	for ( int i = 0; i <= num_max_mismatches; i++) {
		c1 = cnt_k[i];
		for ( int j1 = 0; j1 < nStr; ++j1) {
			for ( int j2 = 0; j2 < nStr; ++j2) {
				K[j1 + j2 * nStr] += w[i] * Ksfinal[(c1 + j1) + j2*nStr];
			}
		}
	}
	/*Normalize kernel values and write into a file*/
	printf("Writing kernel matrix to %s...\n", opfilename);
	FILE *kernelfile;
	FILE *labelfile;
	kernelfile = fopen(opfilename, "w");
	labelfile = fopen(filename_label, "w");
	float kernel_val;

	for (int i = 0; i < nStr; ++i) {
		for (int j = 0; j < nStr; ++j) {
			kernel_val = K[i + j*nStr] / sqrt(K[i + i*nStr] * K[j + j*nStr]);
			if (kernel_val < 1e-50) kernel_val = 0.0;
			fprintf(kernelfile, "%d:%e ", j + 1, kernel_val);
		}
		fprintf(labelfile, "%d ", seqLabels[i]);
		fprintf(labelfile, "\n");
		fprintf(kernelfile, "\n");
	}
	printf("Done\n");
		
	fclose(kernelfile);
	fclose(labelfile);
	free(cnt_k);
	free(Ksfinal);
	free(seqLabels);
	free(K);
	free(nchoosekmat);
	free(feat);
	free(elems);
	return 0;
}
