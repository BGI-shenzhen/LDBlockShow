#ifndef LDDecay_H_
#define LDDecay_H_

#include "HeadIN.h"
#include "FileDeal.h"
#include "Calculate.h"
#include "PairWiseRR.h"
#include "PairWiseDD.h"
#include "PairWiseBoth.h"
#include "GetFig.h"
//#include "SpineBlockRR.h"

using namespace std;

void  LDdecaySNP_help()
{
	cout <<""
		"\n"
		"\tUsage: LDBlockShow  -InVCF  <in.vcf.gz>  -OutPut <outPrefix>  -Region  chr1:10000-20000\n"
		"\n"
		"\t\t-InVCF        <str>      Input SNP VCF Format\n"
		"\t\t-OutPut       <str>      OutPut File of LD Blocks\n"
		"\t\t-Region       <str>      In One Region to show LD info svg Figure\n"
		//"\t\t-NoShow                Only quickly calculate the Spine blcoks\n"
		//"\t\t-Bed        <str>      In muti Region.bed file to show LD\n"
		"\n"
		"\n"
		"\t\t-SeleVar      <int>      Select statistic for deal. 1: D' 2: R^2 3/4:Both [1]\n"
		"\t\t-SubPop       <str>      SubGroup Sample File List[ALLsample]\n"
		"\t\t-BlockType    <int>      method to detect Block [beta] [1]\n"
		"\t\t                         1. Block by PLINK (Gabriel method)\n"
		"\t\t                         2. Solid Spine of LD RR/D' 3. Blockcut with self-defined RR/D'\n"
		"\t\t                         4. FixBlock by input blocks files 5. No Block\n"
		//"\n"
		"\n"
		"\t\t-InGWAS       <str>      InPut GWAS Pvalue File(chr site Pvalue)\n"
		"\t\t-InGFF        <str>      InPut GFF3 file to show Gene CDS and name\n"
		"\n"
		"\t\t-BlockCut     <float>    'Strong LD' cutoff and ratio for BlockType3[0.85:0.90]\n"
		"\t\t-FixBlock     <str>      Input fixed block region\n"
		"\t\t-MerMinSNPNum <int>      merger color grids when SNPnumber over N[50]\n"
		//"\t\t-Ratio       <float>    Ratio of Gabriel Blocks for 'Strong LD' pairwise SNP [0.90]\n"
		//"\t\t-SBMinLen    <float>    Min Length for samall Blocks [10]\n"

		"\t\t\n"
		"\t\t-help                    Show more Parameters and help [hewm2008 v1.36]\n"
		"\n";
}

void More_HelpLDDecay()
{
	cout<<""
		"\n"
		"\t\t More Help document please see the Manual.pdf file\n"
		"\t\t Para [-i] is show for [-InVCF], Para [-o] is show for [-OutPut], Para [-r] is show for [-Region],\n"
		"\n"
		"\n"
		"\t\t-InGenotype    <str>      InPut SNP Genotype Format\n"
		"\t\t-InPlink       <str>      InPut Plink [bed+bim+fam] or [ped+map] file prefix\n"
		"\n"
		"\t\t-MAF           <float>    Min minor allele frequency filter [0.05]\n"
		"\t\t-Het           <float>    Max ratio of het allele filter [0.90]\n"
		"\t\t-Miss          <float>    Max ratio of miss allele filter [0.25]\n"
		"\n"
		"\t\t-TagSNPCut    <float>     'Strong LD' cutoff for TagSNP [0.80]\n"
		"\n"
		"\t\t-OutPng                   convert svg 2 png file\n"
		"\t\t-OutPdf                   convert svg 2 png file\n"
		//		"\t\t-Ratio       <float>    Ratio of Gabriel Blocks for 'Strong LD' pairwise SNP [0.90]\n"
		//		"\t\t-SBMinLen    <float>    Min Length for samall Blocks [10]\n"
		"\n";
}

int LDdecay_help01(int argc, char **argv , In3str1v * paraFA04, Para_18 * para_18)
{
	if (argc <2 ) {LDdecaySNP_help();return 0;}
	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "InVCF" ||  flag  == "i")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr1=argv[i];
		}
		else if (flag  == "SubPop"  ||  flag  == "s")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->SubPop=argv[i];
		}
		else if (flag  == "InGenotype")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag  == "InPlink")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InPlink=argv[i];
		}
		else if (flag  == "InGWAS")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-InGWAS\t"+argv[i];
		}
		else if (flag  == "Cutline")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-Cutline\t"+argv[i];
		}
		else if (flag  == "NoShowLDist")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-NoShowLDist\t"+argv[i];
		}
		else if (flag  == "SpeSNPName")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-SpeSNPName\t"+argv[i];
		}
		else if (flag  == "MerMinSNPNum")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-MerMinSNPNum\t"+argv[i];
		}
		else if (flag  == "InGFF")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-InGFF\t"+argv[i];
		}
		else if (flag  == "FixBlock")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->FixBlock=argv[i];
		}
		else if (flag  ==  "OutPut" ||  flag  == "o")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
		}
		else if (flag  == "Het")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			para_18->Het=atof(argv[i]);
		}
		else if (flag  == "OutPng")
		{
			(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-OutPng";
		}
		else if (flag  == "ShowGWASSpeSNP")
		{
			(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-ShowGWASSpeSNP";
		}
		else if (flag  == "OutPdf")
		{
			(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-OutPdf";
		}
		else if (flag  == "NoGeneName")
		{
			(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-NoGeneName";
		}
		else if (flag  == "ShowNum")
		{
			(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-ShowNum";
		}
		else if (flag  == "TopSite")
		{
			if(i + 1== argc)
			{
				(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-TopSite";
			}
			else if (argv[i+1][0] == '-' )
			{
			(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-TopSite";
			}
			else
			{
				i++;
				(paraFA04->InGWAS)=(paraFA04->InGWAS)+"\t-TopSite\t"+argv[i];
			}
		}
		else if (flag == "MAF")
		{
			if(i + 1== argc) {LogLackArg(flag);return 0;}
			i++;
			para_18->MAF=atof(argv[i]);
		}
		else if (flag == "SBMinLen")
		{
			if(i + 1== argc) {LogLackArg(flag);return 0;}
			i++;
			paraFA04->SBMinLen=atoi(argv[i]);
		}
		else if (flag == "BlockCut")
		{
			if(i + 1== argc) {LogLackArg(flag);return 0;}
			i++;
			string  tmpStr=argv[i];
			if (tmpStr.find(':') ==string::npos)
			{					
				paraFA04->MinRR=atof(argv[i]);
			}
			else
			{
				//Stat.substr(Stat.rfind('/') ==string::npos ? Stat.length() : Stat.rfind('/') + 1);
				string tmpRatioRR=tmpStr.substr(tmpStr.find(':')+1);
				string tmpMinRR=tmpStr.substr(0,tmpStr.find(':'));
				paraFA04->MinRR  = atof(tmpMinRR.c_str());
				paraFA04->RatioRR=atof(tmpRatioRR.c_str());
				if  ( (paraFA04->RatioRR)>1  || (paraFA04->RatioRR)<0.2 )
				{
					cerr <<"\t\tRatio value  should be [0.2,1]\n";
					return 0;
				}
			}
			if  ( (paraFA04->MinRR)>1  || (paraFA04->MinRR)<0.2 )
			{
				cerr <<"\t\t-BlockCut  should be [0.2,1]\n";
				return 0;
			}
		}
		else if (flag == "Ratio")
		{
			if(i + 1== argc) {LogLackArg(flag);return 0;}
			i++;
			paraFA04->RatioRR=atof(argv[i]);
			if  ( (paraFA04->RatioRR)>1  || (paraFA04->RatioRR)<0 )
			{
				cerr <<"\t\t-Ratio  should be (0,1]\n";
				return 0;
			}
		}
		else if (flag == "TagSNPCut")
		{
			if(i + 1== argc) {LogLackArg(flag);return 0;}
			i++;
			paraFA04->TagRR=atof(argv[i]);
			if  ( (paraFA04->TagRR)>1  || (paraFA04->TagRR)<0 )
			{
				cerr <<"\t\t-TagSNPCut  should be (0,1]\n";
				return 0;
			}
		}
		else if (flag == "Miss")
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			para_18->Miss=atof(argv[i]);
		}
		else if (flag == "Region" ||   flag == "r" )
		{
			if(i + 1 == argc) {LogLackArg(flag);return 0;}
			i++;
			paraFA04->Region=argv[i];
			if ((paraFA04->Region).find(':') ==string::npos)
			{
				cerr<<"\tPara [-Region] should be [chr:Start:End],such [chr1:5000:8000]"<<endl;
				return 0;
			}
		}		
		/*
		   else if (flag == "Bed" )
		   {
		   if(i + 1 == argc) {LogLackArg(flag);return 0;}
		   i++;
		   paraFA04->Bed=argv[i];
		   }
		   */		
		   else if (flag == "BlockType")
		   {
			   if(i + 1 == argc) {LogLackArg(flag);return 0;}
			   i++;
			   paraFA04->Method=atoi(argv[i]);
		   }
		   else if (flag == "SeleVar")
		   {
			   if(i + 1 == argc) {LogLackArg(flag);return 0;}
			   i++;
			   paraFA04->RR_D=atoi(argv[i]);
		   }
		   else if (flag == "help"  ||  flag == "h")
		   {
			   More_HelpLDDecay();return 0;
		   }
		   else if (flag == "NoShow")
		   {
			   paraFA04->NoShow=false;
		   }	
		   else
		   {
			   cerr << "UnKnow argument -"<<flag<<endl;
			   return 0;
		   }
	}

	if  ( (paraFA04->InStr2).empty()  )
	{
		cerr<< "lack argument for the must [-OutPut]"<<endl;
		return 0;
	}

	if ( (paraFA04->InStr1).empty() &&  (paraFA04->InStr3).empty() &&  (paraFA04->InPlink).empty() )
	{
		cerr<< "lack argument for InPut File..."<<endl;
		return 0;
	}

	if  ( ((paraFA04->Region).empty())  &&  (paraFA04->NoShow)  )
	{
		cerr<< "\t\tPara [-Region]  should give out like chr1:1-100000\n";
		return 0;
	}

	if ( ((paraFA04->Method)==1) &&   ((para_18->MAF)<0.01) )
	{
		(para_18->MAF)=0.01;
		cerr<<"\t\t -MAF for [-BlockType 1] should >=0.05 ; but now  -MAF reset  0.01 Now"<<endl;
	}

	if ( (paraFA04->Method) >5  || (paraFA04->Method) <0 )
	{
		cerr<<"\t\t-BlockType   should be  [1/2/3/4/5]"<<endl;
		return 0;
	}
	if ( (paraFA04->RR_D) > 4 || (paraFA04->RR_D) < 0 )
	{
		cerr<<"\t\t-SeleVar   should be  [1/2/3/4]"<<endl;
		cerr<<"\t\t1 for D' ; 2 for RR ; 3 for Both RR/D' ; 4 for Both D'/RR\n"<<endl;
		return 0;
	}

	string Stat=(paraFA04->InStr2);
	string ext =Stat.substr(Stat.rfind('.') ==string::npos ? Stat.length() : Stat.rfind('.') + 1);

	if ( ext == "gz" )
	{
		(paraFA04->InStr2)=(paraFA04->InStr2).substr(0,(paraFA04->InStr2).length()-3);
	}

	Stat=(paraFA04->InStr2);
	ext =Stat.substr(Stat.rfind('/') ==string::npos ? Stat.length() : Stat.rfind('/') + 1);

	if (ext != "svg")
	{
		ext =Stat.substr(Stat.rfind('.') ==string::npos ? Stat.length() : Stat.rfind('.') + 1);
		if (ext == "svg")
		{
			(paraFA04->InStr2)=(paraFA04->InStr2).substr(0,(paraFA04->InStr2).length()-4);
		}
	}

	Stat=argv[0];
	GetBinPath ( Stat  ,paraFA04->BinDir);

	return 1 ;
}


int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	Para_18 * para_18 = new Para_18 ;
	if ((LDdecay_help01(argc, argv, paraFA04, para_18)==0))
	{
		delete paraFA04 ;
		delete para_18 ;
		return 1;
	}

	(paraFA04->InInt)=(paraFA04->InInt)*1000;

	//	char buf[1024*1024]; 	setbuf(stdout, buf);
	//*///////////////////////////Test  Out File is OK ///////////////////////////////*//
	string Stat=(paraFA04->InStr2);
	Stat=(paraFA04->InStr2)+".blocks.gz";	
	ogzstream OUTTest ((Stat).c_str());
	if((!OUTTest.good()))
	{
		cerr << "open OUT File error: "<<Stat<<endl;
		delete para_18;
		delete paraFA04 ; return 1;
	}
	OUTTest.close();



	/////////////////////////////////////////////////////////

	map <string,map <llong, vector <BaseType> > >  SNPList ;
	int Flag_for_pro=0;

	///         deal with plink format input ///
	if (!(paraFA04->InPlink).empty())
	{
		if  (PreDealPLink( paraFA04,  para_18 )==0)  /////*   change  ped+map  or bed+bim+fam   2 VCF File /////
		{
			delete para_18;
			delete paraFA04;
			return 1;
		}
	}


	/////*               VCF   IN Deal //////////////////////*////
	if (!(paraFA04->InStr1).empty())
	{
		bool PhaseVCF=false ;
		if  (PreDealVCF( paraFA04,  para_18 ,PhaseVCF )==0)  /////*   PreDeal  VCF  Deal */////
		{
			delete para_18;
			delete paraFA04;
			return 1;
		}

		if (PhaseVCF)
		{
			Read_VCF_IN_Phase( paraFA04, para_18 , SNPList, Flag_for_pro);
		}
		else
		{
			Read_VCF_IN( paraFA04, para_18 , SNPList, Flag_for_pro);
		}
		/*///		
		  if  (PhaseVCF)
		  {		
		  if ((paraFA04->SubPop).empty())
		  {
		  Read_VCF_IN_Phase( paraFA04, para_18 , SNPList, Flag_for_pro);
		  }
		  else
		  {
		  Read_SubPopVCF_IN_Phase( paraFA04, para_18 , SNPList, Flag_for_pro);
		  }
		  }
		  else
		  {
		  if ((paraFA04->SubPop).empty())
		  {
		  Read_VCF_IN( paraFA04, para_18 , SNPList, Flag_for_pro);
		  }
		  else
		  {
		  Read_SubPopVCF_IN( paraFA04, para_18 , SNPList, Flag_for_pro);
		  }
		  }
		  */////
	}

	/////*               Genotype IN Deal //////////////////////*////

	if (!(paraFA04->InStr3).empty())
	{
		if  ((paraFA04->SubPop).empty())
		{
			Read_Genotype_IN(paraFA04, para_18 ,SNPList,Flag_for_pro);
		}
		else
		{
			Read_SubPopGenotype_IN(paraFA04, para_18 ,  SNPList , Flag_for_pro );
		}
	}


	//*///////////////////////////PairWise Compare//////////////////////////////////*//
	if ( (paraFA04->RR_D) ==1 )
	{
		Region_PairWise_D(paraFA04, para_18 ,SNPList, Flag_for_pro);
	}
	else if ( (paraFA04->RR_D) ==2 )
	{
		Region_PairWise_RR(paraFA04, para_18 ,SNPList, Flag_for_pro);
	}
	else if ( (paraFA04->RR_D) ==3 )
	{
		Region_PairWise_BothA(paraFA04, para_18 ,SNPList, Flag_for_pro);
	}
	else
	{
		Region_PairWise_BothB(paraFA04, para_18 ,SNPList, Flag_for_pro);
	}

	if (!(paraFA04->FixBlock).empty())
	{
		AppFixBlock( paraFA04, para_18 );
	}
	//cout<<"Start Draw...\n";
	GetFig_Fun( paraFA04,  para_18 ) ;


	if (!(paraFA04->InStr1).empty())
	{
		Stat="rm -rf "+(paraFA04->InStr1);
		std::system(Stat.c_str()) ;
	}

	delete para_18 ;
	delete paraFA04 ;
	return 0;

}


#endif // LDDecay_H_ //
///////// swimming in the sky and flying in the sea ////////////


