
////////////////////////swimming in the sea & flying in the sky //////////////////


/*
 * DataClass.h
 *
 *  Created on: 2011-11-21
 *      Author: hewm@genomics.org.cn
 */

#ifndef DataClass_H_
#define DataClass_H_


using namespace std;


class In3str1v {
	public:
		string InStr1;
		string InStr2;
		string InStr3;
		string SubPop;
		string InGWAS;
		string Region;
		string BinDir;
		string FixBlock;
		string InPlink;
		int  RR_D ;
		int  InInt ;
		bool NoShow ;
		bool TF2 ;
		int  Method ;
		double TagRR ;
		double MinRR ;
		double RatioRR ;
		int SBMinLen ;
		In3str1v()
		{
			InStr1="";
			InStr2="";
			InStr3="";
			SubPop="";
			InGWAS="";
			InPlink="";
			RR_D=1;
			NoShow=true;
			TF2=true;
			Method=1;
			InInt=300;
			Region="";
			FixBlock="";
			BinDir="";
			TagRR=0.80;
			MinRR=0.90;
			SBMinLen=10;
			RatioRR=0.90;
		}
};


class Para_18 {
	public:
		string input ;
		string output ;
		double Het ;
		double Miss ;
		double MAF ;
		int Cut3base ;
		Para_18()
		{
			input="";
			output="";
			Het=0.90 ;
			Miss=0.25 ;
			Cut3base=0;
			MAF=0.05;
		}
};


class  StarRsult
{
	public:
		int Count;
		double sumRR;
		double sumD;
		StarRsult()
		{
			Count=0;
			sumRR=0.0;
			sumD=0.0;
		}
};


class PairInfoV1
{
	public:
		double D;
		double RR;
		PairInfoV1()
		{
			D=-2;
			RR=0.0;
		}
};


class PairInfoV2
{
	public:
		double D;
		double RR;
		double LOD;
		PairInfoV2()
		{
			D=-2;
			RR=0.0;
			LOD=0.0;
		}
};

class PairInfoV3
{
	public:
		double D;
		double RR;
		double LOD;
		short int low_i ;
		short int high_i ;
		PairInfoV3()
		{
			D=-2;
			RR=0.0;
			LOD=0.0;
			low_i=0;
			high_i=0;
		}
};



class statementVar
{
	public:
		double LN10;
		int Asize ;
		
		double ALL_count ;
		map <string,int>  :: iterator it;
		double D_A;
		double Cal_A ;
		double Cal_B ;
		double D_max;


		double XpA1_pA2;
		double XpA1_pB2;
		double XpB1_pA2;
		double XpB1_pB2;


		double known[5];
		double probHaps[4];
		double lsurface[101];
		double cut5off ;
		unsigned short int DDE[3][3];
		double pA1, pB1, pA2, pB2, loglike1, loglike0;

		double tmpAA, tmpAB, tmpBA, tmpBB, dpr;// tmp2AA, tmp2AB, tmp2BA, tmp2BB;
		int i;
		short int low_i ;
		short int high_i;
		double total_prob;
		double  sum_prob ;
		//		double tmp;//g,h,m,tmp,r;

		statementVar()
		{
			LN10=log(10.0);
//			AA = 0 ; 			AB = 1 ;			BA = 2 ;			BB = 3 ;
//			int AA ;		int AB ;		int BA ;		int BB ;
			total_prob= 0.0;
			sum_prob=0.0;
			low_i = 0;
			high_i = 0;

		}


};


struct BaseType
{
	unsigned short int Value:2 ;
};


#endif /* DataClass_H_ */

//////////////// swimming in the sky and flying in the sea ////////////////
