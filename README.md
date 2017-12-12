# GaKCo-SVM

GaKCo is a fast and naturally parallelizable algorithm for gapped k-mer based string kernel calculation. GaKCo uses associative arrays to calculate the co-occurrence of substrings using cumulative counting. The algorithm easily scales to large dictionary sizes and high numbers of mismatches.
### Reference Paper
Link: [GaKCo: a Fast GApped k-mer string Kernel using COunting](https://arxiv.org/abs/1704.07468)
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
An executable file named `GaKCo` should now be located in your `GaKCo-SVM/bin` directory. If this is not the case, please refer to the trouble shooting section below. 
## Tutorial
### Running GaKCo
GaKCo takes seven arguments: sequence file (fasta), dictionary file, labels file, g, k, kernel file name, multithreading

        Usage: ./GaKCo <sequenceFile> <dictionaryFile> <labelsFile> <g> <k> <kernelFile> <multithreadingEnabled>
        
        Arguments:
            sequenceFile: set of training and testing examples in FASTA format
            dictionaryFile: file containing the alphabet of characters that appear in the sequences (simple text file)
            labelsFile: file to place labels from the examples (simple text file)
            g: parameter for the size of substrings with gaps (or g-mer size)
            k: parameter for the size of non-gapped substrings within the g-mers
            kernelFile: name of the karnel file that will be computed by GaKCo
            multithreading: 1 to enable multithreading, 0 to disable it
Example:
```
    $ ./GaKCo sequences.fasta data/protein.dictionary.txt proteinLabels.txt 7 5 proteinKernel.txt 1
```
This will use the input files to compute a kernel matrix that can be inputted to an SVM classifier.            

#### Limits
GaKCo can currently process up to 50,000 sequences.

### Running GaKCo with the RUN.sh script
You can run GaKCo using the `RUN.sh` script:
```
    $ bash RUN.sh <trainingFile> <testingFile> <dictionaryFile>
```
Where the training and testing files are FASTA files.
Example:
```
    $ bash RUN.sh data/1.1.test.train.fasta data/1.1.train.fasta data/protein.dictionary.txt
```
Alternatively, you can simply use:
```
    $ bash RUN.sh
```
This option uses default hard-coded file names. You can change these default values by opening RUN.sh and changing the file names (located on lines 8, 9, 10)

#### Results
Results (including the computed `kernel.txt` file are placed in your `GaKCo-SVM/results` directory.

## Trouble Shooting
### For Mac Users
Compiling GaKCo requires the OpenMP library, which is not always included with clang (we compile with g++, which usually references clang on Macs). To successfully install GaKCo it may be necessary to do the following:
First check if your machine has Homebrew installed:
```
    $ brew --version
```
If it does not return a version of Homebrew, enter:
```
    $ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```
Once you have verified that Homebrew has been installed, enter:
```
    $ brew install gcc
```
This will create a second GCC installation (that includes OpenMP) located in your /usr/local/bin directory. Open `Makefile` in `GaKCo/src` and edit the CXX line (line 4) to say the following:
```
CXX = g++-7
```
That is, just add "-7" after "g++".
Now save the change and enter:
```
    $ make clean
    $ make all install
```
An executable named "GaKCo" should now be located in your GaKCo/bin directory.

Finally, it may be necessary to install the Command Line Tools package (e.g., if you get an error saying something like "stdio.h not found"):
```
    $ xcode-select --install
```


