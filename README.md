# LDBlockShow
###  1 Introduction
<b>LDBlockShow</b> is a fast and effective tool to generate linkage disequilibrium (LD) heatmap from VCF files. It is more time and memory saving than other current tools. LDBlockShow can generate the plots of LD heatmap and interested statistics or annotation results simultaneously. In addition, it also supports subgroup analysis.

The <i><b>[LDBlockShow article](https://doi.org/10.1093/bib/bbaa227)</b></i> has been published in <b> [briefings in bioinformatics]( https://doi.org/10.1093/bib/bbaa227)</b>, please cited this article if possible

 
PMID: [33126247](https://pubmed.ncbi.nlm.nih.gov/33126247)   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            DOI:[10.1093/bib/bbaa227]( https://doi.org/10.1093/bib/bbaa227)

###  2 Download and Install
------------

<b> 2.1. linux/MaxOS&nbsp;&nbsp;&nbsp;   [Download](https://github.com/BGI-shenzhen/LDBlockShow/archive/v1.36.tar.gz)</b>
  
  </br> <b>2.2 Pre-install</b>
  </br> LDblockshow is for Linux/Unix/macOS only. Before installing, please make sure the following pre-requirements are ready to use.
  </br> 1) g++   : g++ with [--std=c++11](https://gcc.gnu.org/) > 4.8+ is recommended
  </br> 2) zlib  : [zlib](https://zlib.net/) > 1.2.3 is recommended
  </br> 3) Perl  : The [SVG.pm](https://metacpan.org/release/SVG) in Perl should be installed. LDBlockShow uses this module to plot figures. We have provided a built-in SVG module in the package.


</br> <b>2.3 Install</b>
</br> Users can install it with the following options:
</br> Option 1: 
<pre>
        git clone https://github.com/BGI-shenzhen/LDBlockShow.git
        chmod 755 configure  ;  ./configure;
        make;
        mv LDBlockShow  bin/;    #     [rm *.o]
</pre>

</br>**Note:** If failed to link,try to <b>re-install</b> the libraries [**_zlib_**](https://zlib.net/)

</br>Option 2:
<pre>
        tar -zxvf  LDBlockShowXXX.tar.gz
        cd LDBlockShowXXX;  cd src;
        sh  make.sh                         ## Linux :  [ make ; make clean ]
        ../bin/LDBlockShow
</pre>
**Note:** For <b>macOS</b>, if plink doesn't work, please re-download [**mac _plinks_**](https://www.cog-genomics.org/plink2/) and put it into the directory [LDBlockShowXXX/bin]
</br>**Note:** If failed to link,try to <b>re-install</b> the libraries [**_zlib_**](https://zlib.net/)

###  3 Parameter description
------------
</br><b>3.1 LDBlockShow</b>
</br><b>3.1.1 Main parameter</b>

```php
 ./bin/LDBlockShow

        Usage: LDBlockShow  -InVCF  <in.vcf.gz>  -OutPut <outPrefix>  -Region  chr1:10000-20000

                -InVCF          <str>     Input SNP VCF Format
                -OutPut         <str>     OutPut File of LD Blocks
                -Region         <str>     In One Region to show LD info svg Figure


                -SeleVar        <int>     Select statistic for deal. 1: D' 2: R^2 3/4: Both [1]
                -SubPop         <str>     SubGroup Sample File List[ALLsample]
                -BlockType      <int>     method to detect Block [beta] [1]
                                           1. Block by PLINK (Gabriel method)
                                           2. Solid Spine of LD RR/D' 3. Blockcut with self-defined RR/D'
                                           4. FixBlock by input blocks files 5. No Block

                -InGWAS         <str>     InPut GWAS Pvalue File(chr site Pvalue)
                -InGFF          <str>     InPut GFF3 file to show Gene CDS and name

                -BlockCut       <float>   'Strong LD' cutoff and ratio for BlockType3[0.85:0.90]
                -FixBlock       <str>      Input fixed block region
                -MerMinSNPNum   <int>      merger color grids when SNPnumber over N[50]

                -help                      Show more Parameters and help [hewm2008 v1.34]

```
</br> Details for above parameters:
<pre>
-InVCF          The input file in VCF format
-OutPut         The output file directory and output file name prefix (e.g., /path/pop1)
-Region         The defined region to show the LD heatmap (format: chr:start:end)

-SeleVar        The LD measurement (1: D' 2: R^2 3/4: Both R^2 and D'), the default is 1.
-SubPop         A sample list for subgroup analysis
-BlockType      The definition of blocks. The default 1 is called by PLINK1 to generate the block defined by Gabriel et al.2. Solid spine of LD3 is also supported [2]. Users can also define their own cutoff of r2 and D' for blocks [3] combined with the option of "-BlockCut" or supply their own block region definition [4] combined with the option of "-FixBlock". 5 can be used as input if users prefer to not show the block region.

-InGWAS         The statistics file (e.g., association statistics, but other values such as Tajima's D can also be accepted) for generate plot together with the LD plot. File formatted as: [chr position Pvalue]
-InGFF          Input GFF3 format file for genomic region annotation

-BlockCut       For block type 3, the defined cutoff for strong LD, and the ratio of strong LD SNP in one block. Default is 0.85:0.9. That's, if the user chose D' in the -SeleVar option, then in one block, the ratio of SNP pairs with' over 0.85 is 0.9.
-FixBlock       For block type 4, users can use this option to supply a self-defined block region. The file contains three columns, including chromosome, block region start position, and block region end position.
-MerMinSNPNum   The minimum SNP number to merge color grids with the same color. Default is 50. Details please see Fig 1 in this manual.

-help           Show more parameters
</pre>

</br><b>3.1.2 Other parameters</b>
```php
./bin/LDBlockShow -h
		 
		 Para [-i] is show for [-InVCF], Para [-o] is show for [-OutPut], Para [-r] is show for [-Region]

		-InGenotype    <str>     InPut SNP Genotype Format
		-InPlink       <str>     InPut Plink [bed+bim+fam] or [ped+map] file prefix
		
		-MAF           <float>   Min minor allele frequency filter [0.05]
		-Het           <float>   Max ratio of het allele filter [0.90]
		-Miss          <float>   Max ratio of miss allele filter [0.25]

		-TagSNPCut     <float>  'Strong LD' cutoff for TagSNP [0.80]
		-OutPng                  convert svg 2 png file
		-OutPdf                  convert svg 2 pdf file

```
</br> Details for above parameters:
<pre>
-InGenotype   Input file in genotype format.
-InPlink      The prefix of input file in PLINK format.
-MAF          Filter SNPs with low minor allele frequency (default <=0.05)
-Het          Filter SNPs with high heterozygosis ratio (default >=0.9)
-Miss         Filter SNPs with high missing rate (default >=0.25)

-TagSNPCut    The LD cutoff for selecting tag SNPs. Default is 0.8.
-OutPng       Convert the SVG file to PNG file
-OutPdf       Convert SVG file to Pdf file.
</pre>
<b>Note:</b> If users failed to open small SVG files, please use the "-Outpdf" option to use the PDF file. For large SVG files, "-OutPng" can be used to get a relatively small figure file.

</br><b>3.2 ShowLDSVG</b>
</br>This program is designed for users to optimize the figure (e.g., change colors) generated by LDBlockShow.
</br><b>3.2.1 Brief parameters</b>
```php
./bin/ShowLDSVG
		Options
		 -InPreFix    <s> :  InPut Region LD Result Frefix
		 -OutPut      <s> :  OutPut svg file result

		 -help            :  Show more help with more parameter
```
<pre>
-InPreFix     The prefix of input file (i.e., the output file of LDBlockShow)
-OutPut       The out file (svg, png and pdf format plot files)
-help         More parameters in detail
</pre>

</br><b>3.2.2 Detail parameters</b>
```php
./bin/ShowLDSVG   -h
		 -InGWAS      <s>     : InPut GWAS Pvalue File(chr site Pvalue)
		 -NoLogP              : Do not get the log Pvalue
		 -Cutline     <f>     : show the cut off line of Pvalue
		 -TopSite     <n>     : InPut the Special Site as the peak site(chr:pos)
		 -PointSize   <n>     : set the GWAS point size number
		 -SpeSNPName  <s>     : In File for Special SNP Name(chr site Name)
		 -ShowGWASSpeSNP      : show Special SNP Name in GWAS plot with [-SpeSNPName]

		 -InGFF       <s>     : InPut GFF3 file to show Gene CDS and name
		 -NoGeneName          : No show Gene name,only show stuct
		 -crGene      <s>     : InColor for Gene Stuct [CDS:Intron:UTR:Intergenic]
		                        default: [#e7298a:lightblue:#7570b3:#a6cee3]

		 -crBegin     <s>     : In Start Color RGB [255,255,255]
		 -crMiddle    <s>     : In Middle Color RGB [240,235,75]
		 -crEnd       <s>     : In End Color RGB [255,0,0]
		 -NumGradien  <s>     : In Number of gradien of color
		 -crTagSNP    <s>     : Color for TagSNP [31,120,180]

		 -CrGrid      <s>     : the color of grid edge [white]
		 -WidthGrid   <s>     : the edge-width of gird [1]
		 -NoGrid              : No Show the gird edge

		 -ShowNum             : Show the R^2/D' in the heatmap
		 -NoShowLDist   <n>   : NoShow long physical distance pairwise[1000000]
		 -MerMinSNPNum  <s>   : merge color grids when SNPnumber over N[50]

		 -OutPng              : convert svg 2 png file
		 -OutPdf              : convert svg 2 pdf file
		 -ResizeH             : resize image height; Width be resize in ratio[4096]
```
</br> Details for above parameters:
<pre>
-InGWAS         The statistics file (e.g., association statistics, but other values such as Tajima's D can also be accepted) for generate plot together with the LD plot. File formatted as: [chr position Pvalue]
-NoLogP         By default, the P value from the -InGWAS file will be -log10 transformed, with this option, the P value will not be transformed.
-Cutline        The significance cutline of the -InGWAS file.
-TopSite        Users can use this option to assign one interested SNP in the GWAS plot (Default is the most significant SNP, can be changed with chr:pos).
-PointSize      Users can use this option (any number over 0) to set the point size.

-InGFF          The GFF file for genomic region annotation. By default, the gene name will be shown in the plot;
-NoGeneName     Gene name will not be shown in the plot with this option.
-SpeSNPName     With this option, users can input a file to indicate the names for interested SNPs, these names will be shown in the heatmap.
-ShowGWASSpeSNP Users can use this option together with the file assigned by '-SpeSNPName' to show the names of interested SNPs in GWAS plot.
-crGene         Define the colors of different genomic regions. By default, CDS, intron, UTR and intergenic regions will be shown in #e7298a, light blue, #7570b3, and #a6cee3, respectively.


Parameters to optimize the color of the heatmap:
-crBegin        Color for no LD (R^2/D'=0) default: white
-crMiddle       Color for R^2/D'=0.5, default: yellow
-crEnd          Color for complete LD (R^2/D'=1), default: red
-NumGradien     The number of gradients from crBegin to crEnd
-crTagSNP       Color for the tag SNP.

Parameters to optimize the grids in the heatmap:
-CrGrid         Border color of the grids, default: white
-WidthGrid      The width of the border, default = 1
-NoGrid         No border

-ShowNum        Show the LD measurement value in the grids (not recommended when SNP number is over 50).
-NoShowLDist    When the distance between SNPs over this number, their pairwise LD will not be showed in the figure. Default is 10,000,000.
-MerMinSNPNum   When number of SNPs over the default 50, ShowLDSVG will merge adjacent same color grids. User can change this number to any integer numbers.
-OutPng         Convert the SVG file to PNG file
-OutPdf         Convert SVG file to Pdf file.
-ResizeH        Set the height of the image (default 4096), which can be used to adjust the resolution for PNG file. The width will be adjusted automatically.
</pre>

</br><b>Note:</b> When SNP number is large (e.g., over 100), the output SVG file might be very large. ShowLDSVG will merge adjacent same color grids. With smaller number of gradients (set by -NumGradien), the figure will be compressed to be smaller. -MerMinSNPNum is used to set the minimum number of SNPs, that's, if there is more SNPs than this number (default 50), the output SVG will be compressed.

</br><b>3.3 Output files</b>
<pre>
out.site.gz: Remained SNPs after filtering [chr site]
out.blocks.gz: Block file [chr start end block_length SNP_number SNPs]
out.TriangleV.gz: Region Pairwise R^2/D'
out.svg: Output plot in SVG format
out.png: Output plot in png format
out.pdf: Output plot in pdf format
</pre>


###  4 Example
------------

</br>See more detailed usage in the&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <b>[Chinese Documentation](https://github.com/BGI-shenzhen/LDBlockShow/blob/master/LDBlockShow_Manual_Chinese.pdf)</b>
</br>See more detailed usage in the&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <b>[English Documentation](https://github.com/BGI-shenzhen/LDBlockShow/blob/master/LDBlockShow_Manual_English.pdf)</b>
</br>See the example directory and  Manual.pdf for more detail.


* Example 1)  show Figure with Defaut LD Blocks

<pre>
#../../bin/LDBlockShow   -InVCF Test.vcf.gz   -OutPut  out   -Region  chr11:24100000:24200000  -OutPng -SeleVar 1
../../bin/LDBlockShow   -InVCF Test.vcf.gz   -OutPut  out   -Region  chr11:24100000:24200000  -OutPng -SeleVar 2
# [-SeleVar 1] is D',[-SeleVar 2] is RR ,[-SeleVar 3] are RR and D',[-SeleVar 4] are D' and RR # the default is D'
</pre>

![out.png](https://github.com/BGI-shenzhen/LDBlockShow/blob/master/example/Example1/out.png)


* Example 2)   Output LDHeatMap combined with GWAS statistics
```
#../../bin/LDBlockShow  -InVCF ../Example1/Test.vcf.gz -OutPut out -Region chr11:24100000:24200000 -InGWAS gwas.pvalue  -OutPng
../../bin/LDBlockShow   -InVCF ../Example1/Test.vcf.gz -OutPut out -Region chr11:24100000:24200000 -InGWAS gwas.pvalue  -OutPng  -SeleVar 4
## you can run ShowLDSVG with more parameters to optimize the plot  ##
# ../../bin/ShowLDSVG    -InPreFix       out     -OutPut out     -InGWAS gwas.pvalue -Cutline  7  -ShowNum -PointSize 3 -OutPng
```

![out.png](https://github.com/BGI-shenzhen/LDBlockShow/blob/master/example/Example2/out.png)


* Example 3) show Figure with genomic annotation
```
#../../bin/LDBlockShow   -InVCF   ../Example1/Test.vcf.gz    -OutPut  out -InGWAS gwas.pvalue  -InGFF  In.gff   -Region  chr11:24100000:24200000 -OutPng -SeleVar 1   
../../bin/LDBlockShow   -InVCF   ../Example1/Test.vcf.gz    -OutPut  out -InGWAS gwas.pvalue  -InGFF  In.gff   -Region  chr11:24100000:24200000 -OutPng -SeleVar 2
## you can run ShowLDSVG with more parameters to optimize the plot  ##
#../../bin/ShowLDSVG    -InPreFix out -OutPut out.svg -InGWAS gwas.pvalue  -Cutline  7  -InGFF  In.gff  -crGene yellow:lightblue:pink:orange -showNum -OutPng
#../../bin/ShowLDSVG    -InPreFix       out     -OutPut out.svg         -InGFF  In.gff
#../../bin/ShowLDSVG    -InPreFix out -OutPut out.svg -InGWAS gwas.pvalue  -Cutline  7  -InGFF  In.gff  -crGene yellow:lightblue:pink:orange -showNum -OutPng  -SpeSNPName Spe.snp  -ShowGWASSpeSNP
```

![out.png](https://github.com/BGI-shenzhen/LDBlockShow/blob/master/example/Example3/out.png)


* Example 4) show Figure (heatmap+Annotation+GWAS similar to LocusZoom)
```
../../bin/LDBlockShow   -InVCF   ../Example1/Test.vcf.gz    -OutPut  out -InGWAS ../Example3/gwas.pvalue  -InGFF  ../Example3/In.gff   -Region  chr11:24100000:24200000 -OutPng -SeleVar 4 -TopSite
#  [-SeleVar 3]: GWAS with RR ,heatmap with D'.  [-SeleVar 4]: GWAS with D' ,heatmap with RR.
## you can run ShowLDSVG with more parameters to optimize the plot with para [-TopSite]  ##
../../bin/ShowLDSVG    -InPreFix out -OutPut out.svg -InGWAS ../Example3/gwas.pvalue  -Cutline  7  -InGFF  ../Example3/In.gff  -crGene yellow:lightblue:pink:orange -showNum -OutPng  -SpeSNPName ../Example3/Spe.snp  -ShowGWASSpeSNP  -TopSite
#../../bin/ShowLDSVG    -InPreFix out -OutPut out.svg -InGWAS ../Example3/gwas.pvalue  -Cutline  7  -InGFF ../Example3/In.gff  -crGene yellow:lightblue:pink:orange -showNum -OutPng  -SpeSNPName ../Example3/Spe.snp  -ShowGWASSpeSNP  -TopSite  chr11:24142660

```

![out.png](https://github.com/BGI-shenzhen/LDBlockShow/blob/master/example/Example4/out.png)


###  5 Advantages
To evaluate the performance of LDBlockShow, we used test VCF files to generate the LD heatmap by using LDBlockShow, Haploview[4], LDheatmap[5] and gpart. The calculated r<sup>2</sup> and D’values of LDBlockShow is the same with other tools. As shown in Figure below, LDBlockShow is more time and memory saving than other tools.
</br></br>
![Advantages](https://github.com/BGI-shenzhen/LDBlockShow/blob/master/example/Fig/out2.png)
</br> The above figure shows the comparison of computing cost for LDBlockShow, Ldheatmap, Haploview and gpart. CPU time (A) and memory cost (B) for different methods are shown with a fixed SNP number of 100 and sample size ranging from 2,000 to 60,000. CPU time (C) and memory cost (D) for different methods are shown with a fixed sample size of 2,000 and SNP number ranged from 100 to 1,200. When testing datasets in A-D, both LDBlockShow and gpart finished the analyses within reasonable time and memory. We further tested their performance when handling large dataset. CPU time (E) and memory cost (F) for these two methods are shown with a fixed sample size of 100,000 and SNP number ranged from 300 to 2,500. Computation is performed with one thread of an Intel Xeon CPU E5-2630 v4.
</br></br> As shown in Table below, LDBlockShow can generate the plots of LD heatmap and interested statistics or annotation results simultaneously. In addition, LDBlockShow also supports subgroup analysis.
</br>
|          <b>Performance</b>              |    <b>LDBlockShow</b>   |   <b>Haploview</b>    |  <b>LDheatmap</b>     |      <b>gpart</b>     |
|:-----------------------------------------|:-----------------------:|:---------------------:|:---------------------:|:---------------------:|
|<i><b>Input</b></i>                       |                         |                       |                       |                       |
|Compressed VCF file                       |            &radic;      |        &times;        |         &times;       |         &times;       |
|Uncompressed VCF file                     |            &radic;      |        &times;        |         &times;       |         &radic;       |
|Support subgroup analysis                 |            &radic;      |        &times;        |         &times;       |         &times;       |
|<i><b>Output</b></i>                      |                         |                       |                       |                       |
|Visualize additional statistics           |            &radic;      |        &times;        |         &times;       |         &times;       |
|Visualize genomic annotation              |            &radic;      |        &times;        |         &times;       |         &radic;       |
|Compressed SVG                            |            &radic;      |        &times;        |         &times;       |         &times;       |
|PNG file                                  |            &radic;      |        &radic;        |         &times;       |         &radic;       |
|Block region                              |            &radic;      |        &radic;        |         &times;       |         &radic;       |
|LD measurement                            |     R<sup>2</sup>/D'    |   R<sup>2</sup>/D'    |   R<sup>2</sup>       |    R<sup>2</sup>/D'   | 


###  6 An example image generated by LDBlockShow.
------------

![heatmap](https://github.com/BGI-shenzhen/LDBlockShow/blob/master/example/Fig/out1.png)

###  7 Discussing
------------
- [:email:](https://github.com/BGI-shenzhen/LDBlockShow) hewm2008@gmail.com / hewm2008@qq.com
- join the<b><i> QQ Group : 125293663</b></i>
- The <i><b>[LDBlockShow article](https://doi.org/10.1093/bib/bbaa227)</b></i> has been published in <b> [briefings in bioinformatics]( https://doi.org/10.1093/bib/bbaa227)</b>, please cited this article if possible
- PMID: [33126247](https://pubmed.ncbi.nlm.nih.gov/33126247)   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;            DOI:[10.1093/bib/bbaa227]( https://doi.org/10.1093/bib/bbaa227)

######################swimming in the sky and flying in the sea #############################

