#ifndef ReadDataIn_H_
#define ReadDataIn_H_


#include "./FilterGenotype.h"


using namespace std;

int GetBestBase (string Genotype, map <string ,char > SNP_Allele , map <string,string >  & SNP_back_Allele ,  vector <BaseType> &  Base_list)
{
	vector<string>  Base1 ;
	split(Genotype, Base1," \t");
	map <char,int > Count ;
	int Asize=Base1.size();
	for (int ii=0 ;ii<Asize ; ii++)
	{
		string A_tmp=SNP_back_Allele[Base1[ii]];
		Count[A_tmp[0]]++;
		Count[A_tmp[1]]++;
	}
	char best_base='N';
	char sed_base='N';
	int Max=0;
	int SeD=0;
	map <char,int>  :: iterator it=Count.begin();
	for ( ; it!=Count.end(); it++ )
	{
		if ( (it->first ) == 'N' )
		{
			continue ;
		}
		else if ((it->second)  > Max )
		{
			SeD=Max;
			sed_base=best_base;
			Max=(it->second);
			best_base=it->first;
		}
		else if ( (it->second)  >= SeD )
		{
			SeD=(it->second);
			sed_base=it->first;
		}
	}

	map <string ,string > Allele2double;
	Allele2double["-"]="NN";
	Allele2double["N"]="NN";
	Allele2double["n"]="NN";
	string tmp="";
	string A_base=best_base+tmp;
	string B_base=sed_base+tmp;
	string Het_base=A_base+B_base;
	string C_base=SNP_Allele[Het_base]+tmp;

	Allele2double[A_base]=(A_base+A_base);
	Allele2double[B_base]=(B_base+B_base);
	Allele2double[C_base]=Het_base;
	BaseType TempType ;
	//	TempType.Value=0;		Base_list.push_back(TempType);
	//	TempType.Value=1;		Base_list.push_back(TempType);

	for (int ii=0 ;ii<Asize ; ii++)
	{
		string A_tmp=Allele2double[Base1[ii]];
		if (A_tmp=="") 		{	A_tmp=Het_base;	}

		if (A_tmp[0] == best_base )
		{
			TempType.Value=0;
		}
		else if (A_tmp[0] == sed_base )
		{
			TempType.Value=1;
		}
		else
		{
			TempType.Value=2; 
		}
		Base_list.push_back(TempType);

		if (A_tmp[1] == best_base  )
		{
			TempType.Value=0;
		}
		else if (A_tmp[1] == sed_base )
		{
			TempType.Value=1;
		}
		else
		{
			TempType.Value=2; 
		}
		Base_list.push_back(TempType);
	}


	return 1;
}


/*
   int Read_SubPopVCF_IN(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > > &  SNPList ,int & Flag_for_pro )
   {
   igzstream SampleList ((paraFA04->SubPop).c_str(),ifstream::in);
   if (SampleList.fail())
   {
   cerr << "open Sub Group IN File error: "<<(paraFA04->SubPop)<<endl;
   return  0;
   }

   map <string ,int >  SubVetor;
   map <string ,int >  :: iterator it;

   while(!SampleList.eof())
   {
   string  line ;
   getline(SampleList,line);
   if (line.length()<=0 || line[0] == '#' )  { continue  ; }
   vector<string> inf;
   split(line,inf," \t");
   int A=inf.size();
   for(int ii=0 ; ii<A; ii++)
   {
   it=SubVetor.find(inf[ii]);
   if (it==SubVetor.end())
   {
   SubVetor.insert(map <string ,int> ::value_type(inf[ii],0));
   }
   }
   }
   SampleList.close();

   igzstream VCFIN ((paraFA04->InStr1).c_str(),ifstream::in);
   if (VCFIN.fail())
   {
   cerr << "open VCF File IN File error: "<<(paraFA04->InStr1)<<endl;
   return  0;
   }

   vector <int> SampleSite; 
   vector<string> Vsample ;

   while(!VCFIN.eof())
   {
   string  line ;
   getline(VCFIN,line);
   if (line.length()<=0 )  { continue  ; }
   else if ( line[0] == '#' && line[1] == '#' )  { continue  ; }
   else if ( line[0] == '#' && line[1] != '#' )
   {
   Vsample.clear();
   split(line,Vsample," \t");
   if  ( Vsample[0]  != "#CHROM")
   {
   continue  ;
   }
   int A=Vsample.size();

   for (int ii=9 ; ii< A ; ii++)
   {
   it=SubVetor.find(Vsample[ii]);
   if (it!=SubVetor.end())
   {
   SampleSite.push_back(ii);
   (it->second)++;
   }
   }
   break ;
   }
   else if ( line[0] != '#' && line[1] != '#' )
   {
cerr<<"wrong Line : "<<line<<endl;
cerr<<"VCF Header something wrong, can't find sample info before site info"<<endl;
cerr<<"VCF Header sample info Flag : [  #CHROM  ] "<<endl;
return  0;
break;
}
}

int NumberSubGroup=SampleSite.size();
cout<<"the Number of subPop samples[found in VCF] is "<<NumberSubGroup<<endl;
if (NumberSubGroup<3)
{
	cerr<<"sub Group Population szie is too small, SubGroup sample size: "<<NumberSubGroup<<endl;
	return  0;
}


for(it=SubVetor.begin(); it!=SubVetor.end() ;  it++)
{
	if ((it->second)==0)
	{
		cerr<<"warning : Sample [ "<<(it->first)<<" ] can't be found in the VCF Header\n";
	}
	else if ((it->second)>1)
	{
		cerr<<"warning : Sample [ "<<(it->first)<<" ] can be found [Repeat] in the VCF Header\n";
	}
}


int Asample=Vsample.size();

int BadSite=0;
int BadIndelSite=0;

vector<string> Alt ;
vector<string> Btmp  ;
int  Base_len ;
llong Site ;
int Het_count=0;
int Miss_count=0;
char ABase ;
char BBase ;

int BaseConut=0;
char best_base='N';
char sed_base='N';
int Max=0;
int SeD=0;
map <char,int>  :: iterator itSS ;
int ERA=0;
vector<char>  genotype ;
vector<BaseType>  genotypeVE ;


while(!VCFIN.eof())
{
	string  line ;
	getline(VCFIN,line);
	if (line.length()<=0  )  { continue  ; }
	istringstream isoneLine (line,istringstream::in);
	for (int iik=0 ; iik<Asample ; iik++)
	{
		isoneLine>>Vsample[iik];
	}

	Base_len=Vsample[3].length();
	Alt.clear();
	split(Vsample[4],Alt,",");
	for (int ii=0 ; ii<Alt.size() ;ii++)
	{
		if (Alt[ii].length()>Base_len)
		{
			Base_len=Alt[ii].length();
		}
	}

	if (Base_len>1)
	{
		BadIndelSite++;
		continue ;
	}

	map <char,int > Count ;
	Het_count=0;
	Miss_count=0;

	for (int kk=0 ; kk< NumberSubGroup ; kk++)
	{
		Btmp.clear();
		split(Vsample[SampleSite[kk]], Btmp,":");
		string Genotype=Btmp[0];

		if (  Genotype[0] == '.' )
		{
			Miss_count++ ;
		}
		else
		{
			if  (Genotype[0] != Genotype[2] )
			{
				Het_count++;
			}
			Count[Genotype[0]]++;
			Count[Genotype[2]]++;
		}
	}


	if ( ( (Miss_count*1.0/NumberSubGroup)  >(para_18->Miss)  )  )
	{
		continue ;
	}

	if ( ( (Het_count*1.0/NumberSubGroup)  >(para_18->Het) )  )
	{
		continue ;
	}

	BaseConut=0;
	best_base='N';
	sed_base='N';
	Max=0;
	SeD=0;		

	for ( itSS=Count.begin();  itSS!=Count.end(); itSS++ )
	{
		if ( (itSS->first ) == 'N' )
		{
			continue ;
		}
		else if ((itSS->second)  > Max )
		{
			SeD=Max;
			sed_base=best_base;
			Max=(itSS->second);
			best_base=itSS->first;
		}
		else if ( (itSS->second)  >= SeD )
		{
			SeD=(itSS->second);
			sed_base=itSS->first;
		}
		BaseConut++;
	}
	if (BaseConut==1 || BaseConut >2 )
	{
		BadSite++;
		continue ;
	}

	//if ( (  (1-(Max*0.5/NumberSubGroup))  < (para_18->MAF) )  )
	if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
	{
		continue ;
	}

	genotypeVE.clear();
	genotype.clear();
	//		vector<BaseType>  genotypeVE ;

	BaseType  TypeA;

	//		vector<char>  genotype ;
	for (int kk=0 ; kk< NumberSubGroup ; kk++)
	{
		Btmp.clear();
		split(Vsample[SampleSite[kk]], Btmp,":");
		string Genotype=Btmp[0];
		ABase=Genotype[0];
		if (  ABase == '.' )
		{
			genotype.push_back('N');
			genotype.push_back('N');
		}
		else
		{
			BBase=Genotype[2];
			if  (ABase != BBase)
			{
				genotype.push_back(best_base);   // best and sed base to phase
				genotype.push_back(sed_base);
			}
			else
			{
				genotype.push_back(ABase);
				genotype.push_back(BBase);
			}
		}
	}


	ERA=genotype.size();
	for (int hh=0 ; hh<ERA ;hh++)
	{
		if (genotype[hh] == best_base )
		{
			TypeA.Value=0;
		}
		else if (genotype[hh] == sed_base )
		{
			TypeA.Value=1;
		}
		else
		{
			TypeA.Value=2; 
		}
		genotypeVE.push_back(TypeA);	
	}

	istringstream isone (Vsample[1],istringstream::in);
	isone>> Site ;


	map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(Vsample[0]);

	if (itSNP == SNPList.end())
	{
		map <llong, vector <BaseType> > DD;
		DD[Site]=genotypeVE;
		SNPList.insert(map <string,map <llong,vector <BaseType> > > ::value_type(Vsample[0],DD));
		Flag_for_pro++;
	}
	else
	{
		(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE)) ;
		Flag_for_pro++;
	}
}

if(BadIndelSite!=0)
{
	cout<<"warning skip Indel site, there are total skip Indel sites number is : "<<BadIndelSite<<endl;
}
if (BadSite!=0)
{
	cout<<"Warning skip non bi-allelic(Singleton/ThreeMulti allelic) site, and total skip allelic sites number is :"<<BadSite<<endl;
}

return 1;
}
///*////

/*/////////////////
  bool DeleVCFPhase( string  VCFfile )
  {

  igzstream VCFIN ( VCFfile.c_str(),ifstream::in);
  if (VCFIN.fail())
  {
  cerr << "open VCF File IN File error: "<<VCFfile<<endl;
  return  false;
  }

  vector<string> inf ;
  bool  TTFF=false ;
  while(!VCFIN.eof())
  {
  string  line ;
  getline(VCFIN,line);
  if (line.length()<=0 )  { continue  ; }
  else if ( line[0] == '#'  )  { continue  ; }
  else
  {
  split(line,inf," \t");
  if ((inf[9])[1] == '|') 
  {
  TTFF=true;
  }
  break ;
  }
  }
  VCFIN.close();
  if  (TTFF)
  {
  cout <<"#Detected VCF File is phased file with '|', Read VCF in Phase mode"<<endl;
  }
  return TTFF;
  }
///*/////////////





int Read_VCF_IN(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > > &  SNPList ,int & Flag_for_pro )
{

	igzstream VCFIN ((paraFA04->InStr1).c_str(),ifstream::in);

	if (VCFIN.fail())
	{
		cerr << "open VCF File IN File error: "<<(paraFA04->InStr1)<<endl;
		return  0;
	}

	vector<string> inf ;

	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0 )  { continue  ; }
		else if ( line[0] == '#' && line[1] == '#' )  { continue  ; }
		else if ( line[0] == '#' && line[1] != '#' )
		{
			inf.clear();
			split(line,inf," \t");
			if  ( inf[0]  != "#CHROM")
			{
				continue  ;
			}			
			break ;
		}
		else if ( line[0] != '#' && line[1] != '#' )
		{
			cerr<<"wrong Line : "<<line<<endl;
			cerr<<"VCF Header something wrong, can't find sample info before site info"<<endl;
			cerr<<"VCF Header sample info Flag : [  #CHROM  ] "<<endl;
			return  0;
			break;
		}
	}

	int BadSite=0 ;
	int BadIndelSite=0;

	int Asample=inf.size();
	int SampleNum=(Asample-9);


	vector<string> Alt ;
	vector<string> Btmp  ;
	int  Base_len ;
	llong Site ;
	int Het_count=0;
	int Miss_count=0;

	int BaseConut=0;
	char best_base='N';
	char sed_base='N';
	int Max=0;
	int SeD=0;
	map <char,int>  :: iterator it ;
	int ERA=0;
	vector<char>  genotype ;
	vector<BaseType>  genotypeVE ;

	if (paraFA04->TF2)
	{

		while(!VCFIN.eof())
		{
			string  line ;
			getline(VCFIN,line);
			if (line.length()<=0)  { continue  ; }

			istringstream isoneLine (line,istringstream::in);
			for (int iik=0 ; iik<Asample ; iik++)
			{
				isoneLine>>inf[iik];
			}
			Base_len=inf[3].length();
			Alt.clear();
			split(inf[4],Alt,",");

			for (int ii=0 ; ii<Alt.size() ;ii++)
			{
				if (Alt[ii].length()>Base_len)
				{
					Base_len=Alt[ii].length();
				}
			}

			if (Base_len>1)
			{
				BadIndelSite++;
				continue ;
			}

			map <char,int > Count ;
			Het_count=0;
			Miss_count=0;

			for (int jj=9 ; jj< Asample ;jj++ )
			{
				Btmp.clear();
				split(inf[jj], Btmp,":");
				string Genotype=Btmp[0];
				if (  Genotype[0] == '.' )
				{
					Miss_count++ ;
				}
				else
				{
					if  (Genotype[0] != Genotype[2] )
					{
						Het_count++;
					}
					Count[Genotype[0]]++;
					Count[Genotype[2]]++;
				}
			}

			if ( ( (Miss_count*1.0/SampleNum)  >(para_18->Miss)  )  )
			{
				continue ;
			}

			if ( ( (Het_count*1.0/SampleNum)  >(para_18->Het) )  )
			{
				continue ;
			}

			BaseConut=0;
			best_base='N';
			sed_base='N';
			Max=0;
			SeD=0;		

			for ( it=Count.begin(); it!=Count.end(); it++ )
			{
				if ( (it->first ) == 'N' )
				{					
					continue ;
				}
				else if ((it->second)  > Max )
				{
					SeD=Max;
					sed_base=best_base;
					Max=(it->second);
					best_base=it->first;
				}
				else if ( (it->second)  >= SeD )
				{
					SeD=(it->second);
					sed_base=it->first;
				}
				BaseConut++;
			}
			if (BaseConut==1 || BaseConut >2 )
			{
				BadSite++;
				continue ;
			}

			if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
			{
				continue ;
			}

			genotype.clear();

			for (int jj=9 ; jj< Asample ;jj++ )
			{
				Btmp.clear();
				split(inf[jj], Btmp,":");
				string Genotype=Btmp[0];
				if (  Genotype[0] == '.' )
				{
					genotype.push_back('N');
					genotype.push_back('N');
				}
				else
				{
					if  (Genotype[0] != Genotype[2] )	
					{
						genotype.push_back(best_base);
						genotype.push_back(sed_base);
					}	
					else
					{
						genotype.push_back(Genotype[0]);
						genotype.push_back(Genotype[2]);
					}
				}
			}



			genotypeVE.clear();

			BaseType  TypeA;

			ERA=genotype.size();

			for (int hh=0 ; hh<ERA ;hh++)
			{

				if (genotype[hh] == best_base  )
				{
					TypeA.Value=0;
				}
				else if (genotype[hh] == sed_base )
				{
					TypeA.Value=1;
				}
				else
				{
					TypeA.Value=2; 
				}
				genotypeVE.push_back(TypeA);
			}





			istringstream isone (inf[1],istringstream::in);
			isone>> Site ;


			map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(inf[0]);

			if (itSNP == SNPList.end())
			{
				map <llong, vector <BaseType> > DD;
				DD[Site]=genotypeVE;
				SNPList.insert(map <string,map <llong,vector <BaseType> > > ::value_type(inf[0],DD));
				Flag_for_pro++;
			}
			else
			{
				(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE)) ;
				Flag_for_pro++;
			}
		}



	}
	else
	{

		string OUT_VCFTMP=(paraFA04->InStr2)+".vcf.filter.gz";
		ogzstream OUTVCFFF ((OUT_VCFTMP).c_str());

		while(!VCFIN.eof())
		{
			string  line ;
			getline(VCFIN,line);
			if (line.length()<=0 || line[0] == '#' )  { continue  ; }
			llong Site ;
			//inf.clear();			split(line,inf," \t");
			istringstream isoneLine (line,istringstream::in);
			for (int iik=0 ; iik<Asample ; iik++)
			{
				isoneLine>>inf[iik];
			}
			Base_len=inf[3].length();
			Alt.clear();
			split(inf[4],Alt,",");
			for (int ii=0 ; ii<Alt.size() ;ii++)
			{
				if (Alt[ii].length()>Base_len)
				{
					Base_len=Alt[ii].length();
				}
			}

			if (Base_len>1)
			{
				BadIndelSite++;
				continue ;
			}

			map <char,int > Count;
			Het_count=0;
			Miss_count=0;

			for (int jj=9 ; jj< Asample ;jj++)
			{
				Btmp.clear();
				split(inf[jj], Btmp,":");
				string Genotype=Btmp[0];
				if (  Genotype[0]  == '.' )
				{
					Miss_count++ ;
				}
				else
				{
					if  (Genotype[0] != Genotype[2] )
					{
						Het_count++;
					}
					Count[Genotype[0]]++;
					Count[Genotype[2]]++;
				}
			}

			//			int SampleNum=(Asample-9);
			if ( ( (Miss_count*1.0/SampleNum)  >(para_18->Miss)  )  )
			{
				continue ;
			}

			if ( ( (Het_count*1.0/SampleNum)  >(para_18->Het) ) )
			{
				continue ;
			}

			BaseConut=0;
			best_base='N';
			sed_base='N';
			Max=0;
			SeD=0;

			for ( it=Count.begin(); it!=Count.end(); it++ )
			{
				if ( (it->first ) == 'N' )					
				{
					continue ;
				}
				else if ((it->second)  > Max )
				{
					SeD=Max;
					sed_base=best_base;
					Max=(it->second);
					best_base=it->first;
				}
				else if ( (it->second)  >= SeD )
				{
					SeD=(it->second);
					sed_base=it->first;
				}
				BaseConut++;
			}
			if (BaseConut==1 || BaseConut >2 )
			{
				BadSite++;
				continue ;
			}

			//if ( (  (1-(Max*1.0/(SeD+Max)))  < (para_18->MAF) )  )
			if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
			{
				continue ;
			}


			genotype.clear();
			for (int jj=9 ; jj< Asample ;jj++ )
			{
				Btmp.clear();
				split(inf[jj], Btmp,":");
				string Genotype=Btmp[0];
				if (  Genotype[0] == '.' )
				{
					genotype.push_back('N');
					genotype.push_back('N');
				}
				else
				{
					if  (Genotype[0] != Genotype[2] )
					{
						genotype.push_back(best_base);
						genotype.push_back(sed_base);
					}
					else
					{
						genotype.push_back(Genotype[0]);
						genotype.push_back(Genotype[2]);
					}
				}
			}



			//vector<BaseType>  genotypeVE ;
			genotypeVE.clear();

			BaseType  TypeA;

			ERA=genotype.size();
			for (int hh=0 ; hh<ERA ;hh++)
			{
				if (genotype[hh] == best_base  )
				{
					TypeA.Value=0;
				}
				else if (genotype[hh] == sed_base )
				{
					TypeA.Value=1;
				}
				else
				{
					TypeA.Value=2; 
				}
				genotypeVE.push_back(TypeA);
			}


			istringstream isone (inf[1],istringstream::in);
			isone>> Site ;


			map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(inf[0]);
			if (itSNP == SNPList.end())
			{
				map <llong, vector <BaseType> > DD;
				DD[Site]=genotypeVE;
				SNPList.insert(map <string,map <llong,vector <BaseType> > > ::value_type(inf[0],DD));
				Flag_for_pro++;
			}
			else
			{
				(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE)) ;
				Flag_for_pro++;
			}

			OUTVCFFF<<line<<"\n";

		}


		OUTVCFFF.close();

	}
	VCFIN.close();

	if(BadIndelSite!=0)
	{
		cout<<"warning skip Indel site, there are total skip Indel sites number is : "<<BadIndelSite<<endl;
	}
	if (BadSite!=0)
	{
		cout<<"Warning skip non bi-allelic(Singleton/ThreeMulti allelic) site, and total skip allelic sites number is :"<<BadSite<<endl;
	}

	return 1;
}




int Read_Genotype_IN(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > > &  SNPList,int &  Flag_for_pro )
{
	string OUT_TMP=(paraFA04->InStr2)+".genotype.filter.gz";

	char * A = const_cast<char*>((paraFA04->InStr3).c_str());
	char * B = const_cast<char*>((OUT_TMP).c_str());
	stringstream   sstrmC ;
	sstrmC  <<  (para_18->MAF);
	string  C=sstrmC.str();
	char * MAF = const_cast<char*>((C).c_str());

	stringstream   sstrmA ;
	sstrmA  <<  (para_18->Miss);
	string  D=sstrmA.str();
	char * Miss = const_cast<char*>((D).c_str());

	stringstream   sstrmB ;
	sstrmB  <<  (para_18->Het);
	string	 E=sstrmB.str();
	char * Het = const_cast<char*>((E).c_str());

	string  Para_1="FilterGeno";
	char * str_Para_1=const_cast<char*>((Para_1).c_str());

	string	Para_2="-InPut";
	char * str_Para_2=const_cast<char*>((Para_2).c_str());

	string	Para_3="-OutPut";
	char * str_Para_3=const_cast<char*>((Para_3).c_str());

	string	Para_4="-Cut3base";
	char * str_Para_4=const_cast<char*>((Para_4).c_str());

	string	Para_5="-Miss";
	char * str_Para_5=const_cast<char*>((Para_5).c_str());

	string	Para_6="-MAF";
	char * str_Para_6=const_cast<char*>((Para_6).c_str());

	string	Para_7="-Het";
	char * str_Para_7=const_cast<char*>((Para_7).c_str());

	char * TmpFF[12]={str_Para_1, str_Para_2, A , str_Para_3 , B , str_Para_4, str_Para_5, Miss, str_Para_6, MAF ,str_Para_7, Het };
	Filter_genotype_main( 12 ,  TmpFF );

	map <string ,char > SNP_Allele ;
	SNP_Allele["AC"]='M'; SNP_Allele["CA"]='M'; SNP_Allele["GT"]='K'; SNP_Allele["TG"]='K';
	SNP_Allele["CT"]='Y'; SNP_Allele["TC"]='Y'; SNP_Allele["AG"]='R'; SNP_Allele["GA"]='R';
	SNP_Allele["AT"]='W'; SNP_Allele["TA"]='W'; SNP_Allele["CG"]='S'; SNP_Allele["GC"]='S';
	SNP_Allele["AA"]='A'; SNP_Allele["TT"]='T'; SNP_Allele["CC"]='C'; SNP_Allele["GG"]='G';

	map <string,string > SNP_back_Allele ;
	SNP_back_Allele["M"]="AC";SNP_back_Allele["K"]="GT";SNP_back_Allele["Y"]="CT";
	SNP_back_Allele["R"]="AG";SNP_back_Allele["W"]="AT";SNP_back_Allele["S"]="CG";
	SNP_back_Allele["C"]="CC";SNP_back_Allele["G"]="GG";SNP_back_Allele["T"]="TT";
	SNP_back_Allele["A"]="AA";
	SNP_back_Allele["-"]="NN"; SNP_back_Allele["N"]="NN";

	igzstream SNP (OUT_TMP.c_str(),ifstream::in);
	if (SNP.fail())
	{
		cerr << "open SNP File error: "<<OUT_TMP<<endl;
		return  0;
	}


	//////// swimming in the sea & flying in the sky /////////////

	vector<string> inf ;
	while(!SNP.eof())
	{
		string  line ;
		getline(SNP,line);
		if (line.length()<=0 || line[0] == '#' )  { continue  ; }
		llong Site ;
		inf.clear();
		split(line,inf,"\t");

		istringstream isone (inf[1],istringstream::in);
		isone>> Site ;
		vector <BaseType>  genotype ;
		GetBestBase ( inf[2] , SNP_Allele , SNP_back_Allele , genotype ) ;
		map <string,map <llong, vector <BaseType> > >  :: iterator it=SNPList.find(inf[0]);
		if (it == SNPList.end())
		{
			map <llong, vector <BaseType> > DD;
			DD[Site]=genotype;
			SNPList.insert(map <string,map <llong,vector <BaseType> > > ::value_type(inf[0],DD));
			Flag_for_pro++;
		}
		else
		{
			(it->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotype)) ;
			Flag_for_pro++;
		}
	}

	SNP.close();

	if (paraFA04->TF2)
	{
		string  MV="rm -rf  "+OUT_TMP;
		std::system(MV.c_str()) ;
	}


	return 1;
}





int Read_SubPopGenotype_IN(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > > &  SNPList ,int & Flag_for_pro )
{
	igzstream SampleList ((paraFA04->SubPop).c_str(),ifstream::in);
	if (SampleList.fail())
	{
		cerr << "open Sub Group IN File error: "<<(paraFA04->SubPop)<<endl;
		return  0;
	}

	map <string ,int >  SubVetor;
	map <string ,int >  :: iterator it;

	while(!SampleList.eof())
	{
		string  line ;
		getline(SampleList,line);
		if (line.length()<=0 || line[0] == '#' )  { continue  ; }
		vector <string> inf ;
		split(line,inf," \t");
		int A=inf.size();
		for(int ii=0 ; ii<A ; ii++)
		{
			it=SubVetor.find(inf[ii]);
			if (it==SubVetor.end())
			{
				SubVetor.insert(map <string ,int> ::value_type(inf[ii],0));
			}
		}
	}
	SampleList.close();

	igzstream GenotypeIN ((paraFA04->InStr3).c_str(),ifstream::in);
	if (GenotypeIN.fail())
	{
		cerr << "open Genotype File IN File error: "<<(paraFA04->InStr3)<<endl;
		return  0;
	}

	vector <int> SampleSite;

	//vector<string> Vsample ;
	vector <string> inf ;

	while(!GenotypeIN.eof())
	{
		string  line ;
		getline(GenotypeIN,line);
		if (line.length()<=0 )  { continue  ; }
		else if (line[0] == '#' && line[1] == '#' )  { continue  ; }
		else if( line[0] == '#' && line[1] != '#')
		{
			inf.clear();
			split(line,inf," \t");
			if  ( inf[0]  != "#CHROM")  
			{
				continue  ;
			}
			int A=inf.size();
			for(int ii=2 ; ii< A ; ii++)
			{
				it=SubVetor.find(inf[ii]);
				if (it!=SubVetor.end())
				{
					SampleSite.push_back(ii);
					(it->second)++;
				}
			}
			break ;
		}
		else if ( line[0] != '#' && line[1] != '#' )
		{
			cerr<<"wrong Line : "<<line<<endl;
			cerr<<"Genotype Header something wrong, can't find sample info before site info"<<endl;
			cerr<<"iTools   Formtools  VCF2Genotype  -InPut  in.vcf  -OutPut  out.genotype  -WithHeader   -NoRef"<<endl;
			cerr<<"Genotype Header sample info Flag : [  #CHROM  ] "<<endl;
			return  0;
			break;
		}
	}

	int NumberSubGroup=SampleSite.size();
	cout<<"the Number of subPop samples[found in Genotype] is "<<NumberSubGroup<<endl;
	if (NumberSubGroup<3)
	{
		cerr<<"sub Group Population szie is too small, SubGroup sample size: "<<NumberSubGroup<<endl;
		return  0;
	}


	for(it=SubVetor.begin(); it!=SubVetor.end() ;  it++)
	{
		if ((it->second)==0)
		{
			cerr<<"warning : Sample [ "<<(it->first)<<" ] can't be found in the Genotype Header\n";
		}
		else if ((it->second)>1)
		{
			cerr<<"warning : Sample [ "<<(it->first)<<" ] can be found [Repeat] in the Genotype Header\n";
		}
	}


	int BadSite=0;


	map <string,string > SNP_back_Allele ;
	SNP_back_Allele["M"]="AC";SNP_back_Allele["K"]="GT";SNP_back_Allele["Y"]="CT";
	SNP_back_Allele["R"]="AG";SNP_back_Allele["W"]="AT";SNP_back_Allele["S"]="CG";
	SNP_back_Allele["C"]="CC";SNP_back_Allele["G"]="GG";SNP_back_Allele["T"]="TT";
	SNP_back_Allele["A"]="AA";
	SNP_back_Allele["-"]="NN"; SNP_back_Allele["N"]="NN";



	int Asample=inf.size();		


	vector<char>  genotype ;
	vector<BaseType>  genotypeVE ;

	while(!GenotypeIN.eof())
	{
		string  line ;
		getline(GenotypeIN,line);
		if (line.length()<=0  )  { continue  ; }
		llong Site ;
		//	inf.clear();		split(line,inf," \t");
		istringstream isoneLine (line,istringstream::in);
		for (int iik=0 ; iik<Asample ; iik++)
		{
			isoneLine>>inf[iik];
		}


		map <char,int > Count ;
		int Het_count=0;
		int Miss_count=0;



		for (int kk=0 ; kk< NumberSubGroup ; kk++)
		{
			string Genotype=SNP_back_Allele[inf[SampleSite[kk]]];
			char ABase=Genotype[0];
			if (  ABase == 'N' )
			{
				Miss_count++ ;
			}
			else
			{
				char BBase=Genotype[1];
				if  (ABase != BBase )
				{
					Het_count++;
				}
				Count[ABase]++;
				Count[BBase]++;
			}
		}


		if ( ( (Miss_count*1.0/NumberSubGroup)  >(para_18->Miss)  )  )
		{
			continue ;
		}

		if ( ( (Het_count*1.0/NumberSubGroup)  >(para_18->Het) )  )
		{
			continue ;
		}

		int BaseConut=0;
		char best_base='N';
		char sed_base='N';
		int Max=0;
		int SeD=0;		
		map <char,int>  :: iterator it=Count.begin();

		for ( ; it!=Count.end(); it++ )
		{
			if ( (it->first ) == 'N' )
			{
				continue ;
			}
			else if ((it->second)  > Max )
			{
				SeD=Max;
				sed_base=best_base;
				Max=(it->second);
				best_base=it->first;
			}
			else if ( (it->second)  >= SeD )
			{
				SeD=(it->second);
				sed_base=it->first;
			}
			BaseConut++;
		}
		if (BaseConut==1 || BaseConut >2 )
		{
			BadSite++;
			continue ;
		}

		//if ( (  (1-(Max*1.0/(SeD+Max)))  < (para_18->MAF) )  )
		if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
		{
			continue ;
		}


		genotypeVE.clear();
		BaseType  TypeA;

		//		cerr<<best_base<<"\t"<<Max<<"\t"<<sed_base<<"\t"<<SeD<<endl;

		genotype.clear();
		for (int kk=0 ; kk< NumberSubGroup ; kk++)
		{
			string Genotype=SNP_back_Allele[inf[SampleSite[kk]]];
			char ABase=Genotype[0];
			if ( ABase == 'N' )
			{
				genotype.push_back('N');
				genotype.push_back('N');
			}
			else
			{
				char BBase=Genotype[1];
				if  (ABase != BBase )
				{
					genotype.push_back(best_base);
					genotype.push_back(sed_base);
				}
				else
				{
					genotype.push_back(ABase);
					genotype.push_back(BBase);
				}
			}
		}





		int ERA=genotype.size();
		for (int hh=0 ; hh<ERA ;hh++)
		{

			if (genotype[hh] == best_base )
			{
				TypeA.Value=0;
			}
			else if (genotype[hh] == sed_base )
			{
				TypeA.Value=1;
			}
			else
			{
				TypeA.Value=2; 
			}
			genotypeVE.push_back(TypeA);
			//	cerr<<TypeA.Value<<" ";
		}
		//		cerr<<endl;


		istringstream isone (inf[1],istringstream::in);
		isone>> Site ;


		map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(inf[0]);

		if (itSNP == SNPList.end())
		{
			map <llong, vector <BaseType> > DD;
			DD[Site]=genotypeVE;
			SNPList.insert(map <string,map <llong,vector <BaseType> > > ::value_type(inf[0],DD));
			Flag_for_pro++;
		}
		else
		{
			(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE));
			Flag_for_pro++;
		}
	}


	if (BadSite!=0)
	{
		cout<<"Warning skip non bi-allelic(Singleton/ThreeMulti allelic) site, and total skip allelic sites number is :"<<BadSite<<endl;
	}

	return 1;
}







/////////////////////   Phase  VCF /////////////////////

/*////
  int Read_SubPopVCF_IN_Phase(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > > &  SNPList ,int & Flag_for_pro )
  {
  igzstream SampleList ((paraFA04->SubPop).c_str(),ifstream::in);
  if (SampleList.fail())
  {
  cerr << "open Sub Group IN File error: "<<(paraFA04->SubPop)<<endl;
  return  0;
  }

  map <string ,int >  SubVetor;
  map <string ,int >  :: iterator it;

  while(!SampleList.eof())
  {
  string  line ;
  getline(SampleList,line);
  if (line.length()<=0 || line[0] == '#' )  { continue  ; }
  vector<string> inf;
  split(line,inf," \t");
  int A=inf.size();
  for(int ii=0 ; ii<A; ii++)
  {
  it=SubVetor.find(inf[ii]);
  if (it==SubVetor.end())
  {
  SubVetor.insert(map <string ,int> ::value_type(inf[ii],0));
  }
  }
  }
  SampleList.close();

  igzstream VCFIN ((paraFA04->InStr1).c_str(),ifstream::in);
  if (VCFIN.fail())
  {
  cerr << "open VCF File IN File error: "<<(paraFA04->InStr1)<<endl;
  return  0;
  }

  vector <int> SampleSite; 
  vector<string> Vsample ;

  while(!VCFIN.eof())
  {
  string  line ;
  getline(VCFIN,line);
  if (line.length()<=0 )  { continue  ; }
  else if ( line[0] == '#' && line[1] == '#' )  { continue  ; }
  else if ( line[0] == '#' && line[1] != '#' )
  {
  Vsample.clear();
  split(line,Vsample," \t");
  if  ( Vsample[0]  != "#CHROM")
  {
  continue  ;
  }
  int A=Vsample.size();

  for (int ii=9 ; ii< A ; ii++)
  {
  it=SubVetor.find(Vsample[ii]);
  if (it!=SubVetor.end())
  {
  SampleSite.push_back(ii);
  (it->second)++;
  }
  }
  break ;
  }
  else if ( line[0] != '#' && line[1] != '#' )
  {
cerr<<"wrong Line : "<<line<<endl;
cerr<<"VCF Header something wrong, can't find sample info before site info"<<endl;
cerr<<"VCF Header sample info Flag : [  #CHROM  ] "<<endl;
return  0;
break;
}
}

int NumberSubGroup=SampleSite.size();
cout<<"the Number of subPop samples[found in VCF] is "<<NumberSubGroup<<endl;
if (NumberSubGroup<3)
{
	cerr<<"sub Group Population szie is too small, SubGroup sample size: "<<NumberSubGroup<<endl;
	return  0;
}


for(it=SubVetor.begin(); it!=SubVetor.end() ;  it++)
{
	if ((it->second)==0)
	{
		cerr<<"warning : Sample [ "<<(it->first)<<" ] can't be found in the VCF Header\n";
	}
	else if ((it->second)>1)
	{
		cerr<<"warning : Sample [ "<<(it->first)<<" ] can be found [Repeat] in the VCF Header\n";
	}
}


int Asample=Vsample.size();

int BadSite=0;
int BadIndelSite=0;

vector<string> Alt ;
vector<string> Btmp  ;
int  Base_len ;
llong Site ;
int Het_count=0;
int Miss_count=0;
char ABase ;
char BBase ;

int BaseConut=0;
char best_base='N';
char sed_base='N';
int Max=0;
int SeD=0;
map <char,int>  :: iterator itSS ;
int ERA=0;
vector<char>  genotype ;
vector<BaseType>  genotypeVE ;


while(!VCFIN.eof())
{
	string  line ;
	getline(VCFIN,line);
	if (line.length()<=0  )  { continue  ; }
	istringstream isoneLine (line,istringstream::in);
	for (int iik=0 ; iik<Asample ; iik++)
	{
		isoneLine>>Vsample[iik];
	}

	Base_len=Vsample[3].length();
	Alt.clear();
	split(Vsample[4],Alt,",");
	for (int ii=0 ; ii<Alt.size() ;ii++)
	{
		if (Alt[ii].length()>Base_len)
		{
			Base_len=Alt[ii].length();
		}
	}

	if (Base_len>1)
	{
		BadIndelSite++;
		continue ;
	}

	map <char,int > Count ;
	Het_count=0;
	Miss_count=0;

	for (int kk=0 ; kk< NumberSubGroup ; kk++)
	{
		Btmp.clear();
		split(Vsample[SampleSite[kk]], Btmp,":");
		string Genotype=Btmp[0];

		if (  Genotype[0] == '.' )
		{
			Miss_count++ ;
		}
		else
		{
			if  (Genotype[0] != Genotype[2] )
			{
				Het_count++;
			}
			Count[Genotype[0]]++;
			Count[Genotype[2]]++;
		}

	}


	if ( ( (Miss_count*1.0/NumberSubGroup)  >(para_18->Miss)  )  )
	{
		continue ;
	}

	if ( ( (Het_count*1.0/NumberSubGroup)  >(para_18->Het) )  )
	{
		continue ;
	}

	BaseConut=0;
	best_base='N';
	sed_base='N';
	Max=0;
	SeD=0;		

	for ( itSS=Count.begin();  itSS!=Count.end(); itSS++ )
	{
		if ( (itSS->first ) == 'N' )
		{
			continue ;
		}
		else if ((itSS->second)  > Max )
		{
			SeD=Max;
			sed_base=best_base;
			Max=(itSS->second);
			best_base=itSS->first;
		}
		else if ( (itSS->second)  >= SeD )
		{
			SeD=(itSS->second);
			sed_base=itSS->first;
		}
		BaseConut++;
	}

	if (BaseConut==1 || BaseConut >2 )
	{
		BadSite++;
		continue ;
	}

	//if ( (  (1-(Max*0.5/NumberSubGroup))  < (para_18->MAF) )  )
	if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
	{
		continue ;
	}

	genotypeVE.clear();
	genotype.clear();


	BaseType  TypeA;

	for (int kk=0 ; kk< NumberSubGroup ; kk++)
	{
		Btmp.clear();
		split(Vsample[SampleSite[kk]], Btmp,":");
		string Genotype=Btmp[0];
		ABase=Genotype[0];
		if (  ABase == '.' )
		{
			genotype.push_back('N');
			genotype.push_back('N');
		}
		else
		{
			BBase=Genotype[2];
			genotype.push_back(ABase);	// phase VCF
			genotype.push_back(BBase);  // phase VCF
		}
	}


	ERA=genotype.size();
	for (int hh=0 ; hh<ERA ;hh++)
	{
		if (genotype[hh] == best_base )
		{
			TypeA.Value=0;
		}
		else if (genotype[hh] == sed_base )
		{
			TypeA.Value=1;
		}
		else
		{
			TypeA.Value=2;
		}
		genotypeVE.push_back(TypeA);
	}



	istringstream isone (Vsample[1],istringstream::in);
	isone>> Site ;


	map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(Vsample[0]);

	if (itSNP == SNPList.end())
	{
		map <llong, vector <BaseType> > DD;
		DD[Site]=genotypeVE;
		SNPList.insert(map <string,map <llong,vector <BaseType> > > ::value_type(Vsample[0],DD));
		Flag_for_pro++;
	}
	else
	{
		(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE)) ;
		Flag_for_pro++;
	}
}

if(BadIndelSite!=0)
{
	cout<<"warning skip Indel site, there are total skip Indel sites number is : "<<BadIndelSite<<endl;
}
if (BadSite!=0)
{
	cout<<"Warning skip non bi-allelic(Singleton/ThreeMulti allelic) site, and total skip allelic sites number is :"<<BadSite<<endl;
}

return 1;
}

////*///////






int Read_VCF_IN_Phase(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > > &  SNPList ,int & Flag_for_pro )
{
	igzstream VCFIN ((paraFA04->InStr1).c_str(),ifstream::in);

	if (VCFIN.fail())
	{
		cerr << "open VCF File IN File error: "<<(paraFA04->InStr1)<<endl;
		return  0;
	}


	vector<string> inf ;

	while(!VCFIN.eof())
	{
		string  line ;
		getline(VCFIN,line);
		if (line.length()<=0 )  { continue  ; }
		else if ( line[0] == '#' && line[1] == '#' )  { continue  ; }
		else if ( line[0] == '#' && line[1] != '#' )
		{
			inf.clear();
			split(line,inf," \t");
			if  ( inf[0]  != "#CHROM")
			{
				continue  ;
			}			
			break ;
		}
		else if ( line[0] != '#' && line[1] != '#' )
		{
			cerr<<"wrong Line : "<<line<<endl;
			cerr<<"VCF Header something wrong, can't find sample info before site info"<<endl;
			cerr<<"VCF Header sample info Flag : [  #CHROM  ] "<<endl;
			return  0;
			break;
		}
	}

	int BadSite=0 ;
	int BadIndelSite=0;

	int Asample=inf.size();
	int SampleNum=(Asample-9);


	vector<string> Alt ;
	vector<string> Btmp  ;
	int  Base_len ;
	llong Site ;
	int Het_count=0;
	int Miss_count=0;

	int BaseConut=0;
	char best_base='N';
	char sed_base='N';
	int Max=0;
	int SeD=0;
	map <char,int>  :: iterator it ;
	int ERA=0;
	vector<char>  genotype ;
	vector<BaseType>  genotypeVE ;

	if (paraFA04->TF2)
	{

		while(!VCFIN.eof())
		{
			string  line ;
			getline(VCFIN,line);
			if (line.length()<=0)  { continue  ; }

			istringstream isoneLine (line,istringstream::in);
			for (int iik=0 ; iik<Asample ; iik++)
			{
				isoneLine>>inf[iik];
			}
			Base_len=inf[3].length();
			Alt.clear();
			split(inf[4],Alt,",");

			for (int ii=0 ; ii<Alt.size() ;ii++)
			{
				if (Alt[ii].length()>Base_len)
				{
					Base_len=Alt[ii].length();
				}
			}

			if (Base_len>1)
			{
				BadIndelSite++;
				continue ;
			}

			map <char,int > Count ;
			Het_count=0;
			Miss_count=0;

			for (int jj=9 ; jj< Asample ;jj++ )
			{
				Btmp.clear();
				split(inf[jj], Btmp,":");
				string Genotype=Btmp[0];
				if (  Genotype[0] == '.' )
				{
					Miss_count++ ;
				}
				else
				{
					if  (Genotype[0] != Genotype[2] )
					{
						Het_count++;
					}
					Count[Genotype[0]]++;
					Count[Genotype[2]]++;
				}
			}

			if ( ( (Miss_count*1.0/SampleNum)  >(para_18->Miss)  )  )
			{
				continue ;
			}

			if ( ( (Het_count*1.0/SampleNum)  >(para_18->Het) )  )
			{
				continue ;
			}

			BaseConut=0;
			best_base='N';
			sed_base='N';
			Max=0;
			SeD=0;		

			for ( it=Count.begin(); it!=Count.end(); it++ )
			{
				if ( (it->first ) == 'N' )
				{					
					continue ;
				}
				else if ((it->second)  > Max )
				{
					SeD=Max;
					sed_base=best_base;
					Max=(it->second);
					best_base=it->first;
				}
				else if ( (it->second)  >= SeD )
				{
					SeD=(it->second);
					sed_base=it->first;
				}
				BaseConut++;
			}
			if (BaseConut==1 || BaseConut >2 )
			{
				BadSite++;
				continue ;
			}

			if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
			{
				continue ;
			}

			genotype.clear();

			for (int jj=9 ; jj< Asample ;jj++ )
			{
				Btmp.clear();
				split(inf[jj], Btmp,":");
				string Genotype=Btmp[0];
				if (  Genotype[0] == '.' )
				{
					genotype.push_back('N');
					genotype.push_back('N');
				}
				else
				{
					genotype.push_back(Genotype[0]);   // phase VCF
					genotype.push_back(Genotype[2]);   // phase VCF
				}
			}



			genotypeVE.clear();

			BaseType  TypeA;

			ERA=genotype.size();

			for (int hh=0 ; hh<ERA ;hh++)
			{

				if (genotype[hh] == best_base  )
				{
					TypeA.Value=0;
				}
				else if (genotype[hh] == sed_base )
				{
					TypeA.Value=1;
				}
				else
				{
					TypeA.Value=2; 
				}
				genotypeVE.push_back(TypeA);
			}





			istringstream isone (inf[1],istringstream::in);
			isone>> Site ;


			map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(inf[0]);

			if (itSNP == SNPList.end())
			{
				map <llong, vector <BaseType> > DD;
				DD[Site]=genotypeVE;
				SNPList.insert(map <string,map <llong,vector <BaseType> > > ::value_type(inf[0],DD));
				Flag_for_pro++;
			}
			else
			{
				(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE)) ;
				Flag_for_pro++;
			}
		}



	}
	else
	{

		string OUT_VCFTMP=(paraFA04->InStr2)+".vcf.filter.gz";
		ogzstream OUTVCFFF ((OUT_VCFTMP).c_str());

		while(!VCFIN.eof())
		{
			string  line ;
			getline(VCFIN,line);
			if (line.length()<=0 || line[0] == '#' )  { continue  ; }
			llong Site ;
			//inf.clear();			split(line,inf," \t");
			istringstream isoneLine (line,istringstream::in);
			for (int iik=0 ; iik<Asample ; iik++)
			{
				isoneLine>>inf[iik];
			}
			Base_len=inf[3].length();
			Alt.clear();
			split(inf[4],Alt,",");
			for (int ii=0 ; ii<Alt.size() ;ii++)
			{
				if (Alt[ii].length()>Base_len)
				{
					Base_len=Alt[ii].length();
				}
			}

			if (Base_len>1)
			{
				BadIndelSite++;
				continue ;
			}

			map <char,int > Count;
			Het_count=0;
			Miss_count=0;

			for (int jj=9 ; jj< Asample ;jj++)
			{
				Btmp.clear();
				split(inf[jj], Btmp,":");
				string Genotype=Btmp[0];
				if (  Genotype[0]  == '.' )
				{
					Miss_count++ ;
				}
				else
				{
					if  (Genotype[0] != Genotype[2] )
					{
						Het_count++;
					}
					Count[Genotype[0]]++;
					Count[Genotype[2]]++;
				}
			}

			//			int SampleNum=(Asample-9);
			if ( ( (Miss_count*1.0/SampleNum)  >(para_18->Miss)  )  )
			{
				continue ;
			}

			if ( ( (Het_count*1.0/SampleNum)  >(para_18->Het) ) )
			{
				continue ;
			}

			BaseConut=0;
			best_base='N';
			sed_base='N';
			Max=0;
			SeD=0;

			for ( it=Count.begin(); it!=Count.end(); it++ )
			{
				if ( (it->first ) == 'N' )					
				{
					continue ;
				}
				else if ((it->second)  > Max )
				{
					SeD=Max;
					sed_base=best_base;
					Max=(it->second);
					best_base=it->first;
				}
				else if ( (it->second)  >= SeD )
				{
					SeD=(it->second);
					sed_base=it->first;
				}
				BaseConut++;
			}
			if (BaseConut==1 || BaseConut >2 )
			{
				BadSite++;
				continue ;
			}

			//if ( (  (1-(Max*1.0/(SeD+Max)))  < (para_18->MAF) )  )
			if ( (SeD*1.0/(SeD+Max))  < (para_18->MAF) )  
			{
				continue ;
			}


			genotype.clear();
			for (int jj=9 ; jj< Asample ;jj++ )
			{
				Btmp.clear();
				split(inf[jj], Btmp,":");
				string Genotype=Btmp[0];
				if (  Genotype[0] == '.' )
				{
					genotype.push_back('N');
					genotype.push_back('N');
				}
				else
				{
					genotype.push_back(Genotype[0]); // phase VCF
					genotype.push_back(Genotype[2]); // phase VCF
				}
			}



			//vector<BaseType>  genotypeVE ;
			genotypeVE.clear();

			BaseType  TypeA;

			ERA=genotype.size();
			for (int hh=0 ; hh<ERA ;hh++)
			{
				if (genotype[hh] == best_base  )
				{
					TypeA.Value=0;
				}
				else if (genotype[hh] == sed_base )
				{
					TypeA.Value=1;
				}
				else
				{
					TypeA.Value=2; 
				}
				genotypeVE.push_back(TypeA);
			}


			istringstream isone (inf[1],istringstream::in);
			isone>> Site ;


			map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(inf[0]);
			if (itSNP == SNPList.end())
			{
				map <llong, vector <BaseType> > DD;
				DD[Site]=genotypeVE;
				SNPList.insert(map <string,map <llong,vector <BaseType> > > ::value_type(inf[0],DD));
				Flag_for_pro++;
			}
			else
			{
				(itSNP->second).insert(map <llong, vector <BaseType> >  :: value_type(Site,genotypeVE)) ;
				Flag_for_pro++;
			}

			OUTVCFFF<<line<<"\n";

		}


		OUTVCFFF.close();

	}
	VCFIN.close();

	if(BadIndelSite!=0)
	{
		cout<<"warning skip Indel site, there are total skip Indel sites number is : "<<BadIndelSite<<endl;
	}
	if (BadSite!=0)
	{
		cout<<"Warning skip non bi-allelic(Singleton/ThreeMulti allelic) site, and total skip allelic sites number is :"<<BadSite<<endl;
	}

	return 1;
}






/////*        PreDeal  VCF  Deal //////////////////////*//

int PreDealVCF(In3str1v * paraFA04, Para_18 * para_18 , bool  &  TTFF  )
{
	map <string ,int >  SubVetor;
	map <string ,int >  :: iterator it;
	vector <int> SampleSite;
	string cutSample="1-9";
	if (!(paraFA04->SubPop).empty())
	{
		igzstream SampleList ((paraFA04->SubPop).c_str(),ifstream::in);
		if (SampleList.fail())
		{
			cerr << "open Sub Group IN File error: "<<(paraFA04->SubPop)<<endl;
			return  0;
		}

		while(!SampleList.eof())
		{
			string  line ;
			getline(SampleList,line);
			if (line.length()<=0 || line[0] == '#' )  { continue  ; }
			vector<string> inf;
			split(line,inf," \t");
			int A=inf.size();
			for(int ii=0 ; ii<A; ii++)
			{
				it=SubVetor.find(inf[ii]);
				if (it==SubVetor.end())
				{
					SubVetor.insert(map <string ,int> ::value_type(inf[ii],0));
				}
			}
		}
		SampleList.close();
	}

	igzstream VCFINTT ((paraFA04->InStr1).c_str(),ifstream::in);
	string TmpVCFIN=(paraFA04->InStr2)+".region.vcf";
	ofstream OUTAATT (TmpVCFIN.c_str());

	if(!OUTAATT.good())
	{
		cerr << "open OUT File error: "<<TmpVCFIN<<endl;
		return 0;
	}

	if (VCFINTT.fail())
	{
		cerr << "open VCF File IN File error: "<<(paraFA04->InStr1)<<endl;
		return  0;
	}

	vector<string> inf ;

	while(!VCFINTT.eof())
	{
		string  line ;
		getline(VCFINTT,line);
		if (line.length()<=0 )  { continue  ; }
		else if ( line[0] == '#' && line[1] == '#' )  {   OUTAATT<<line<<"\n"; continue  ; }
		else if ( line[0] == '#' && line[1] != '#' )
		{
			inf.clear();			
			split(line,inf," \t");
			if  ( inf[0]  != "#CHROM")
			{
				continue  ;
			}
			if ((paraFA04->SubPop).empty())
			{
				OUTAATT<<line<<"\n";
			}
			else
			{
				OUTAATT<<inf[0];
				for (int ii=1 ; ii< 9 ; ii++)
				{
					OUTAATT<<"\t"<<inf[ii];
				}
				int A=inf.size();
				for (int ii=9 ; ii< A; ii++)
				{
					it=SubVetor.find(inf[ii]);
					if (it!=SubVetor.end())
					{
						SampleSite.push_back(ii);
						(it->second)++;
						OUTAATT<<"\t"<<inf[ii];
						int yy=ii+1;
						cutSample=cutSample+","+Int2Str(yy);
					}
				}
				OUTAATT<<"\n";
			}

			getline(VCFINTT,line);
			inf.clear();
			split(line,inf," \t");
			if ((inf[9])[1] == '|') 
			{
				TTFF=true;
				cout <<"#Detected VCF File is phased file with '|', Read VCF in Phase mode"<<endl;
			}
			break ;
		}
		else if ( line[0] != '#' && line[1] != '#' )
		{
			cerr<<"wrong Line : "<<line<<endl;
			cerr<<"VCF Header something wrong, can't find sample info before site info"<<endl;
			cerr<<"VCF Header sample info Flag : [  #CHROM  ] "<<endl;
			return  0;
		}
	}
	VCFINTT.close();
	OUTAATT.close();


	if (!(paraFA04->SubPop).empty())
	{
		int NumberSubGroup=SampleSite.size();
		cout<<"the Number of subPop samples[found in VCF] is "<<NumberSubGroup<<endl;
		if (NumberSubGroup<3)
		{
			cerr<<"sub Group Population szie is too small, SubGroup sample size: "<<NumberSubGroup<<endl;
			return  0;
		}

		for(it=SubVetor.begin(); it!=SubVetor.end() ;  it++)
		{
			if ((it->second)==0)
			{
				cerr<<"warning : Sample [ "<<(it->first)<<" ] can't be found in the VCF Header\n";
			}
			else if ((it->second)>1)
			{
				cerr<<"warning : Sample [ "<<(it->first)<<" ] can be found [Repeat] in the VCF Header\n";
			}
		}

	}

	inf.clear();
	split((paraFA04->Region),inf,":-");
	string Stat=" cat " + (paraFA04->InStr1);
	if (GzipTrue(paraFA04->InStr1))
	{
		Stat=" zcat "+ (paraFA04->InStr1) ;
	}
	if ((paraFA04->SubPop).empty())
	{
		Stat=Stat+" | awk  \'$1==\""+inf[0]+"\"  && $2>="+inf[1]+"  && $2<="+inf[2]+" \'  |  awk   \'$3=$2\'    OFS=\"\\t\"  >> " +TmpVCFIN ;
	}
	else
	{
		Stat=Stat+" | awk  \'$1==\""+inf[0]+"\"  && $2>="+inf[1]+"  && $2<="+inf[2]+" \'  |  awk   \'$3=$2\'   OFS=\"\\t\" | cut -f "+cutSample+"  >> " +TmpVCFIN ;
	}
	std::system(Stat.c_str()) ;
	paraFA04->InStr1=TmpVCFIN;
	return 1;
	/////////////////////////////////////////////////////////

}


void  AppFixBlock(In3str1v * paraFA04, Para_18 * para_18  )
{
	string tmpStat="gzip -d  "+(paraFA04->InStr2)+".blocks.gz";
	std::system(tmpStat.c_str()) ;
	tmpStat=(paraFA04->InStr2)+".blocks";
	ofstream  OUTB (tmpStat.c_str(), ios::out|ios::app);
	OUTB<<"##User InPut Fixed  Blocks"<<endl;
	igzstream INB ((paraFA04->FixBlock).c_str(),ifstream::in);
	while(!INB.eof())
	{
		string  line ;
		getline(INB,line);
		OUTB<<line<<endl;
	}
	OUTB.close();
	INB.close();
	tmpStat="gzip  "+(paraFA04->InStr2)+".blocks";
	std::system(tmpStat.c_str()) ;
}



int PreDealPLink(In3str1v * paraFA04, Para_18 * para_18  )
{

	string bed=(paraFA04->InPlink)+".bed";
	string fam=(paraFA04->InPlink)+".fam";
	string bim=(paraFA04->InPlink)+".bim";

	string ped=(paraFA04->InPlink)+".ped";
	string map=(paraFA04->InPlink)+".map";

	// ../../bin/plink   --file www  --export vcf --out   yyy  --allow-extra-chr  
	// ../../bin/plink   --bfile www  --export vcf --out   yyyzzz  --allow-extra-chr
	vector<string> RegionBed;
	split((paraFA04->Region),RegionBed,":-");
	string ForRegion=" --chr  "+RegionBed[0]+" --from-bp  "+RegionBed[1] + " --to-bp "+ RegionBed[2];
//--chr  Ghir_D11  --from-bp   24100000     --to-bp  24185529  

	string commod;
	if  ( (access(bed.c_str(), 0) == 0)  &&  (access(fam.c_str(), 0) == 0)  &&  (access(bim.c_str(), 0) == 0)  )
	{
		commod="  --bfile  "+(paraFA04->InPlink)+"  --allow-extra-chr    --export vcf --out  "+(paraFA04->InStr2)+ForRegion ;
	}
	else if  (  (access(map.c_str(), 0) == 0)  &&  (access(ped.c_str(), 0) == 0) )
	{
		commod="  --file  "+(paraFA04->InPlink)+"  --allow-extra-chr    --export vcf --out  "+(paraFA04->InStr2)+ForRegion;
	}
	else
	{
		cerr<<"Can't found These three Files :\n\t"<<bed<<"\n\t"<<fam<<"\n\t"<<bim<<endl;
		cerr<<"And also Can't found These two Files :\n\t"<<ped<<"\n\t"<<map<<endl;
		return 0 ;
	}

	string Plink=(paraFA04->BinDir)+"/plink";

	if ( access(Plink.c_str(), 0) == 0 )
	{
		Plink=Plink+commod +"  > /dev/null ";
		std::system(Plink.c_str());
		Plink="rm -rf "+(paraFA04->InStr2)+".log  " + (paraFA04->InStr2)+".nosex  > /dev/null  "  ;
		std::system(Plink.c_str());
	}
	else
	{
		char   buf[2048]={'\0'};
		string cc="which  plink  2> /dev/null ";
		memset( buf, '\0', sizeof(buf) );
		FILE   *stream ;
		stream=popen(cc.c_str(),"r") ;
		fread( buf, sizeof(char), sizeof(buf), stream);
		Plink=buf;
		Plink=Plink.substr(0,Plink.length()-1);

		if (Plink == "" )
		{
			cout <<"\twarning: can't find the [plink] in your $PATH ; please find it and copy into "<<endl;
			return 0 ;
		}
		else
		{
			Plink=Plink+commod +"  > /dev/null ";
			std::system(Plink.c_str());
			Plink="rm -rf "+(paraFA04->InStr2)+".log  " + (paraFA04->InStr2)+".nosex  > /dev/null  " ;
			std::system(Plink.c_str());
		}
	}

	(paraFA04->InStr1)=(paraFA04->InStr2)+".vcf";
	return 1;
}



////////////

#endif // LDDecay_H_ //
///////// swimming in the sky and flying in the sea ////////////




