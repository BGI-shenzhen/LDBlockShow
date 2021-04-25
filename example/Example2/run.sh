#!/bin/sh
#$ -S /bin/sh
#Version1.0	hewm@genomics.org.cn	2020-02-07
echo Start Time : 
date
#../../bin/LDBlockShow  -InVCF ../Example1/Test.vcf.gz -OutPut out -Region chr11:24100000:24200000 -InGWAS gwas.pvalue  -OutPng
../../bin/LDBlockShow   -InVCF ../Example1/Test.vcf.gz -OutPut out -Region chr11:24100000:24200000 -InGWAS gwas.pvalue  -OutPng  -SeleVar 2
## you can run ShowLDSVG with more parameter to re-plot the svg if you want to make some changes ##
# ../../bin/ShowLDSVG    -InPreFix    out   -OutPut out    -InGWAS gwas.pvalue -Cutline  7  -ShowNum -PointSize 3 -OutPng
echo End Time : 
date
