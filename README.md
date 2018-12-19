# GaKCo-SVM

GaKCo is a fast and naturally parallelizable algorithm for gapped k-mer based string kernel calculation. GaKCo uses associative arrays to calculate the co-occurrence of substrings using cumulative counting. The algorithm easily scales to large dictionary sizes and high numbers of mismatches.
### Reference Paper
Link: [GaKCo: a Fast GApped k-mer string Kernel using COunting](https://arxiv.org/abs/1704.07468)

### Citations

```
@inproceedings{singh_gakco:_2017,
    location = {Cham},
    title = {{GaKCo}: A Fast Gapped k-mer String Kernel Using Counting},
    isbn = {978-3-319-71249-9},
    abstract = {String Kernel ({SK}) techniques, especially those using gapped k-mers as features (gk), have obtained great success in classifying sequences like {DNA}, protein, and text. However, the state-of-the-art gk-{SK} runs extremely slow when we increase the dictionary size (\$\${\textbackslash}{backslashvarSigma} \$\$) or allow more mismatches (M). This is because current gk-{SK} uses a trie-based algorithm to calculate co-occurrence of mismatched substrings resulting in a time cost proportional to \$\$O({\textbackslash}{backslashvarSigma} {\textasciicircum}\{M\})\$\$. We propose a fast algorithm for calculating Gapped k-mer Kernel using Counting ({GaKCo}). {GaKCo} uses associative arrays to calculate the co-occurrence of substrings using cumulative counting. This algorithm is fast, scalable to larger \$\${\textbackslash}{backslashvarSigma} \$\$and M, and naturally parallelizable. We provide a rigorous asymptotic analysis that compares {GaKCo} with the state-of-the-art gk-{SK}. Theoretically, the time cost of {GaKCo} is independent of the \$\${\textbackslash}{backslashvarSigma} {\textasciicircum}\{M\}\$\$term that slows down the trie-based approach. Experimentally, we observe that {GaKCo} achieves the same accuracy as the state-of-the-art and outperforms its speed by factors of 2, 100, and 4, on classifying sequences of {DNA} (5 datasets), protein (12 datasets), and character-based English text (2 datasets). ({GaKCo} is shared as an open source tool at https://github.com/{QData}/{GaKCo}-{SVM}). Code and data related to this chapter are available at: https://doi.org/10.6084/m9.figshare.5434825.},
    pages = {356--373},
    booktitle = {Machine Learning and Knowledge Discovery in Databases},
    publisher = {Springer International Publishing},
    author = {Singh, Ritambhara and Sekhon, Arshdeep and Kowsari, Kamran and Lanchantin, Jack and Wang, Beilun and Qi, Yanjun},
    editor = {Ceci, Michelangelo and Hollmén, Jaakko and Todorovski, Ljupčo and Vens, Celine and Džeroski, Sašo},
    date = {2017}
}
```

### Included Datasets
We perform 19 different classification tasks to evaluate the performance of GaKCo. The tasks pertain to three categories:

        (1) TF binding site prediction (DNA dataset)
        (2) Remote Protein Homology prediction (protein dataset)
        (3) Character-based English text classification (text dataset)
## Installation
Download and extract this repository, then enter:
```
    $ cd src
    $ make all install
```
An executable file named `GaKCo` should now be located in your `GaKCo-SVM/bin` directory. 
## Tutorial
### Running GaKCo
GaKCo takes several arguments:

        Usage: ./GaKCo -g <int> -k <int> -n <int> -p <int> <sequenceFile> <dictionaryFile> <labelsFile> <kernelFile>
        
        Arguments:
            g : length of gapped instance. Constraints: 0 < g < 20
            k : length of k-mer. Constraints: k < g
            n : (optional) maximum number of examples in the data set. Default: 15000 [This value can be increased or decreased according to the memory capacity of the machine]
            p : (optional) parallel. Set to 1 to using multithreading; else set to 0. Default: 1
            sequenceFile : set of sequences (see "Sequence File Format" section below for details on how this file should be formatted)
            dictionaryFile : file containing the alphabet of characters that appear in the sequences (simple text file)
            labelsFile : name of file to place labels from the examples (simple text file). This file is used by GaKCo to create create the kernel matrix
            kernelFile : name of the file to write the kernel that will be computed by GaKCo.
For example:
```
    $ ./GaKCo -g 7 -k 5 -n 15000 -p 1 mySequences.fasta data/protein.dictionary.txt labelsFile.txt computedKernel.txt
```
Or:
```
    $ ./bin/GaKCo -g 7 -k 5 -n 15000 -p 1 mySequences.fasta data/protein.dictionary.txt labelsFile.txt computedKernel.txt
```
#### Output
GaKCo computes the kernel function for all pairs of sequences in the provided sequence file, thus giving the "distance" or "similarity" values needed for an SVM classifier. The resultant kernel matrix will be written to the `kernelFile` argument.

### Running GaKCo with the RUN.sh script
You can run GaKCo using the `RUN.sh` script after completing the "Installation" instructions above. This option is useful if you have separate training and testing files and want to create a training set kernel file and a testing set kernel file.
```
    $ bash RUN.sh <trainingFile.fasta> <testingFile.fasta> <dictionaryFile.txt>
```
Example:
```
    $ bash RUN.sh data/1.1.test.fasta data/1.1.train.fasta data/protein.dictionary.txt
```
Alternatively, you can simply use the following to call the RUN.sh script with default hard-coded values:
```
    $ bash RUN.sh
```
You can change the default values by opening RUN.sh and changing the file names (located on lines 8, 9, 10). You can change the default parameter values by modifying lines 32-35.

#### Output
The script places the computed test and kernel files in the `GaKCo-SVM/results` directory.
### Sequence File Format
The `sequenceFile` argument to GaKCo must be a file in the following format:
```
>label
sequence
>label
sequence
...
```
Where label are class labels (-1 or 0 for negative class and 1 for positive class) and sequence lines are a sequence contained in a single line (sequences running over multiple lines will not be handled correctly). The first line of the file should be a label line. For example:
```
>1
MKTPITEAIAAADNQGRFLSNTELQAVNGRYQRAAASLEAARSL
>1
MLDAFAKVVAQADARGEFLSNTQLDALSKMVSEGNKRLD
>0
FPTIPLSRLADNAWLRADRLNQLAFDTYQEFEEAYIPKEQIHSFWWNPQ
```
