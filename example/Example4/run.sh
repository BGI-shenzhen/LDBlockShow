#!/bin/sh
#$ -S /bin/sh
echo Start Time : 
date
../../bin/LDBlockShow   -InVCF   ../Example1/Test.vcf.gz    -OutPut  out -InGWAS ../Example3/gwas.pvalue  -InGFF  ../Example3/In.gff   -Region  chr11:24100000:24200000 -OutPng -SeleVar 4 -TopSite
#  [-SeleVar 3]: GWAS with RR ,heatmap with D'.  [-SeleVar 4]: GWAS with D' ,heatmap with RR.
##  you can run ShowLDSVG with more parameter to re-plot the svg with para [-TopSite] ##
../../bin/ShowLDSVG    -InPreFix out -OutPut out.svg -InGWAS ../Example3/gwas.pvalue  -Cutline  7  -InGFF  ../Example3/In.gff  -crGene yellow:lightblue:pink:orange -showNum -OutPng  -SpeSNPName ../Example3/Spe.snp  -ShowGWASSpeSNP  -TopSite
#../../bin/ShowLDSVG    -InPreFix out -OutPut out.svg -InGWAS ../Example3/gwas.pvalue  -Cutline  7  -InGFF ../Example3/In.gff  -crGene yellow:lightblue:pink:orange -showNum -OutPng  -SpeSNPName ../Example3/Spe.snp  -ShowGWASSpeSNP  -TopSite  chr11:24142660
date
