#!/bin/sh
#$ -S /bin/sh
echo Start Time : 
date
#../../bin/LDBlockShow   -InVCF   ../Example1/Test.vcf.gz    -OutPut  out -InGWAS gwas.pvalue  -InGFF  In.gff   -Region  chr11:24100000:24200000 -OutPng -SeleVar 1
../../bin/LDBlockShow   -InVCF   ../Example1/Test.vcf.gz    -OutPut  out -InGWAS gwas.pvalue  -InGFF  In.gff   -Region  chr11:24100000:24200000 -OutPng -SeleVar 2
##  you can run ShowLDSVG with more parameter to re-plot the svg if you want to make some changes ##
#../../bin/ShowLDSVG    -InPreFix out -OutPut out.svg -InGWAS gwas.pvalue  -Cutline  7  -InGFF  In.gff  -crGene yellow:lightblue:pink:orange -showNum -OutPng 
#../../bin/ShowLDSVG    -InPreFix       out     -OutPut out.svg         -InGFF  In.gff
# ../../bin/ShowLDSVG    -InPreFix out -OutPut out.svg -InGWAS gwas.pvalue  -Cutline  7  -InGFF  In.gff  -crGene yellow:lightblue:pink:orange -showNum -OutPng  -SpeSNPName Spe.snp  -ShowGWASSpeSNP
date
