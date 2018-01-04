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
GaKCo takes several arguments:

        Usage: ./GaKCo -g <int> -k <int> -n <int> -p <int> <sequenceFile> <dictionaryFile> <labelsFile> <kernelFile>
        
        Arguments:
            g : length of gapped instance. Constraints: 0 < g < 20
            k : length of k-mer. Constraints: k < g
            n : (optional) maximum number of examples in the data set. Default: 15000 [This value can be increased or decreased according                  to the memory capacity of the machine]
            p : (optional) parallel. Set to 1 to using multithreading; else set to 0. Default: 1
            sequenceFile : set of training and testing examples in FASTA format
            dictionaryFile : file containing the alphabet of characters that appear in the sequences (simple text file)
            labelsFile : file to place labels from the examples (simple text file)
            kernelFile : name of the file to write the kernel that will be computed by GaKCo
For example:
```
    $ ./GaKCo -g 7 -k 5 -n 15000 -p 1 mySequences.fasta data/protein.dictionary.txt labelsFile.txt computedKernel.txt
```
This would use the input files to compute a kernel matrix (stored in computedKernel.txt) that can be inputted to an SVM classifier.

### Running GaKCo with the RUN.sh script
You can run GaKCo using the `RUN.sh` script:
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

#### Results
Results (including the computed `kernel.txt` file are placed in your `GaKCo-SVM/results` directory.

## Trouble Shooting
### For Mac Users
(Note: we recommend using a Linux distro as opposed to Mac OSX due to some quirks in the way GCC is usually installed on Macs. But if you strongly prefer using Mac OS X and have run into problems, this section could help.)

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
Now save the change and enter the following while still in the `GaKCo/src` directory:
```
    $ make clean
    $ make all install
```
An executable named "GaKCo" should now be located in your GaKCo/bin directory.

Finally, it may be necessary to install the Command Line Tools package (e.g., if you get an error saying something like "stdio.h not found"):
```
    $ xcode-select --install
```
