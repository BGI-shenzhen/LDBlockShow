#!/usr/bin/env bash
#$ -S /bin/sh
echo Start Time : 
date

echo  Detection system platform  is 
if [ "$(uname)" == "Darwin" ];then
	echo MacOS
	cp  plink_mac  ../bin/plink
	which "gzcat" >/dev/null 2>&1
	if [ $? -eq 0 ];
	then
		echo check gzcat done
		sed  's/ zcat/ gzcat/g' ./FileDeal.h  > ./FileDeal.h.tmp
		mv ./FileDeal.h.tmp  ./FileDeal.h
	else
		echo can not find [gzcat],we used the [zcat] replace,but maybe something wrong
		echo you should install the [gzcat],and then re-install me
	fi
	# Mac OS X 操作系统
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ];then
	echo Linux
	# GNU/Linux操作系统
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ];then
	echo windows
	# Windows NT操作系统
fi

echo g++  -std=c++11 -g -O2 LDBlockShow.cpp		-lm -lc -lz  -L/usr/lib/ -L./include/zlib/	-o	../bin/LDBlockShow
g++ -std=c++11	-g -O2 LDBlockShow.cpp		-lm -lc -lz  -L/usr/lib/ -L./include/zlib/	-o	../bin/LDBlockShow
## for -static 
##  g++	-g -O2 LDBlockShow.cpp		-lm -lc -lz  -L/usr/lib/ -L./include/zlib/	-o	../bin/LDBlockShow  -static
chmod 755 ../bin/* 
echo done see the [ ../bin/LDBlockShow  ]
echo End Time : 
date



