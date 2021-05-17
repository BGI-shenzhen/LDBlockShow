#!/bin/sh
#$ -S /bin/sh
#Version1.0	hewm@genomics.org.cn	2020-01-10
echo Start Time : 
date
g++	-static	-g	-O2	LDBlockShow.cpp	-lz	-L/usr/lib/	-L./include/zlib/	-L	/hwfssz4/BC_PUB/Software/01.Usr/lib/	-o	../bin/LDBlockShow_linux 	-lm	-lc	
chmod 755 ../bin/*
echo End Time : 
date
