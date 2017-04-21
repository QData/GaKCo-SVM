#!/bin/bash
ntrain=`wc -l train.fasta|awk '{print $1/2}'`

ntest=`wc -l test.fasta|awk '{print $1/2}'`
cat train.fasta test.fasta>sequences.fasta

# feed into GaKCo 
#GaKCo <sequencefile> <dictionaryfile> <filename  for labels> <g> <k> <filename for kernel> <set for multithread>   
# g,k (user's choice)
g=7
k=5
GaKCo sequences.fasta dictionary.txt labels.txt $g $k kernel.txt 1

# cut the kernel into train and test
cat kernel.txt|cut -d' ' -f1-$ntrain|head -$ntrain > kernel_train.txt
cat kernel.txt|cut -d' ' -f1-$ntrain|tail -n -$ntest > kernel_test.txt
#cut labels
cat labels.txt|cut -d' ' -f1|head -$ntrain > train.labels.txt
cat labels.txt|cut -d' ' -f1|tail -n $ntest > test.labels.txt


#create empirical feature map (liblinear format : label foloowed by features)
paste -d" " train.labels.txt kernel_train.txt>train.features.txt
paste -d" " test.labels.txt kernel_test.txt>test.features.txt
 #C parameter for liblinear
 C=1
#use liblinear train function (user's choice)
train -c $C -q train.features.txt model.txt
