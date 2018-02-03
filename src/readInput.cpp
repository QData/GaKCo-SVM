// Code Contibution by:
//Ritambhara Singh <rs3zz@virginia.edu>
//Kamran Kowsari <kk7nc@virginia.edu >
//Arshdeep Sekhon <as5cu@virginia.edu >

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>


int * string_replace (char *s, char *d);
int help2();
char * readDict (char *filename, int *na);
int dictsize;

//read input from fasta file

int ** Readinput_(char *filename, char *dictfile, int *Labelout, int *seqLengths, long int *nStr, long int *maxlen, long int *minlen, int *dictionarySize, int max_num_str) {
    int **output;
    char *labelfile, *seqfile;
    char *str, *linetemp, *line, *seq, *trimline, *label;
    int j;
    bool isLabel = true;
    FILE *inpfile;
    char *d;
    int c = 0;
    d = readDict(dictfile, dictionarySize);

    printf("Reading %s\n", filename);
    inpfile = fopen(filename, "r");
    
    if (inpfile) {
        line = (char *) malloc(STRMAXLEN * sizeof(char));
        int row = 0; //counts rows of the output and line number of the sequence file
        output =  (int **) malloc(STRMAXLEN * sizeof(int *));
        bool read = true;
        while (read) {
            printf("=== current row = %d ===\n", row);
            if(!fgets(line, STRMAXLEN, inpfile)) {
                read = false;
            }

            linetemp = (char *) malloc(STRMAXLEN * sizeof(char *));
            strcpy(linetemp, line);
            if (isLabel) {
                label = strtok(linetemp,">");
                if(strcmp(label,"0\n")==0) {
                    strcpy(label,"-1\n");
                }
                Labelout[row]= atoi(label);
                printf("    the label was %d\n", Labelout[row]);
                isLabel = false;
            } else {
                printf("    before trimwhitespace line = %s\n", line);
                trimline = trimwhitespace(line);
                printf("    after trimwhitespace line = %s\n", line);
                strcpy(linetemp, trimline);
                printf("    linetemp = %s\n", linetemp);
                printf("    strlen(linetemp) = %d\n", strlen(linetemp));
                seqLengths[row] = strlen(linetemp);
                printf("    seqLengths[0] = %d\n", seqLengths[0]);
                if (seqLengths[row] > maxlen[0]) {
                    maxlen[0] = seqLengths[row];
                }
                if (seqLengths[row] < minlen[0]) {
                    minlen[0] = seqLengths[row];
                }
                printf("    seqLengths[0] = %d\n", seqLengths[0]);
                output[row] = (int *) malloc(seqLengths[row] * sizeof(int));
                memset(output[row], 0, sizeof(int) * seqLengths[row]);
                strcpy(linetemp, trimline);
                seq = trimline;
                output[row] = string_replace(seq, d);
                row++;
                isLabel = true;
                printf("    seqLengths[0] = %d\n", seqLengths[0]);
            }
            printf("1 seqLengths[0] = %d\n", seqLengths[0]);
            free(linetemp);
            printf("2 seqLengths[0] = %d\n", seqLengths[0]);
        }
        printf("done: seqLengths[0] = %d\n", seqLengths[0]);
        nStr[0] = row;
        fclose(inpfile);
        free(line);
        printf("done: seqLengths[0] = %d\n", seqLengths[0]);
    } else {
        perror(filename);
    }
    printf("\n\n\n");
    for (int x = 0; x < *nStr; x++) {
        printf("seqLengths[%d] = %d\n", x, seqLengths[x]);
    }
    printf("\n\n\n");
    for (int kk = 0; kk < *nStr; kk++) {
        printf("seqLengths[%d] = %d.\n Encoding %d = ", kk, seqLengths[kk], kk);
        for(int jj = 0; jj < seqLengths[kk]; jj++) {
            if(output[kk][jj] > dictsize) {
                output[kk][jj] = 0;
            }
            printf("%d", output[kk][jj]);
        }
        printf("\n");
    }
    return output;
}

// read dictionary to convert into numerical format

char * readDict (char *filename, int *dictionarySize) {
    char *D;
    char *linetemp1, *line1, *next_elem, *trimline;
    int i, j;
    FILE *inpfile;

    inpfile = fopen (filename, "r" );
    D = (char *) malloc(50 * sizeof(char));

    if (inpfile) {
        line1 = (char *) malloc(STRMAXLEN * sizeof(char));
        i = 0;
        while (fgets(line1, STRMAXLEN, inpfile)) {
            linetemp1 = (char *) malloc(STRMAXLEN * sizeof(char));
            trimline = trimwhitespace(line1);
            strcpy(linetemp1, trimline);
            D[i] = linetemp1[0];
            free(linetemp1);
            i++;
        }
        dictsize = i - 1;
        printf("Dictionary size=%d (+1 for uknown character)\n", dictsize + 1);
        fclose(inpfile);
        *dictionarySize = dictsize + 2;
    } else {
        perror(filename);
    }
    return D;
}

//converts g-mers into numerical representation

int * string_replace (char *s, char *d) {
    int i, count, found;
    int *array;
    found = 0;
    array = (int *)malloc(STRMAXLEN*sizeof(int));
    count = 0;
    while(s[count] != '\0') {
        for (i=0; i <= dictsize; i++) {
            if (toupper(s[count]) == d[i]) {
                array[count]=i+1;
                found=1;
            }
        }
        if (found == 0) {
            array[count]=0;
        } else {
            found = 0;
        }
        count++;
    }
    return(array);
}

int help2() {
  printf("Usage: readInput <Input-file> <Labels-file> <Sequence-file> <Dictionary File>\n");
  return 1;
}
