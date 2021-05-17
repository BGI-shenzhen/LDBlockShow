#!/bin/sh
#$ -S /bin/sh
#Version1.0	hewm@genomics.cn	2020-01-10
echo Start Time : 
date
#../../bin/LDBlockShow   -InVCF Test.vcf.gz   -OutPut  out   -Region  chr11:24100000:24200000  -OutPng -SeleVar 1
../../bin/LDBlockShow   -InVCF Test.vcf.gz   -OutPut  out   -Region  chr11:24100000:24200000  -OutPng  -SeleVar 2
# [-SeleVar 1] is D'  , [-SeleVar 2]  is RR . # if you want to use D' for stat ,just defaut
echo End Time : 
date
