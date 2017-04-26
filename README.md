# GaKCo-SVM

Reference Paper: [GaKCo: a Fast GApped k-mer string Kernel using COunting](https://arxiv.org/abs/1704.07468)

GaKCo is a a fast and naturally parallelizable algorithm for gapped k-mer based string kernel calculation. GaKCo uses associative arrays to calculate the co-occurrence of substrings using cumulative counting. This algorithm is scalable to larger dictionary size and more number of mismatches.

**Datasets for GaKCo:** 
We perform 19 different classification tasks to evaluate the performance of GaKCo. These tasks belong to the discussed three categories: (1) TF binding site prediction (DNA dataset), (2) Remote Protein Homology prediction (protein dataset), and (3) Character-based English text classification (text dataset).

The train and test sets have been provided in "data/" folder.

Compiling GaKCo (with openMP) : 
```
g++ -c GaKCo.cpp -o GaKCo -fopenmp
```
To get kernel output: 
```
#GaKCo <sequencefile> <dictionaryfile> <filename  for labels> <g> <k> <filename for kernel> <set for multithread>   
```
Bash script to run end-to-end kernel calculation:
```
processing.sh
```
