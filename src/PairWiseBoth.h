#ifndef PairWiseBoth_H_
#define PairWiseBoth_H_

using namespace std;



int Region_PairWise_BothA(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > >  & SNPList, int & Flag_for_pro)
{
	string StatDD  =(paraFA04->InStr2)+".TriangleV.gz";
	string StatRR  =(paraFA04->InStr2)+".TriangleB.gz";
	string StatSite=(paraFA04->InStr2)+".site.gz";
	string StatBlocks=(paraFA04->InStr2)+".blocks.gz";
	ogzstream OUTDD ((StatDD).c_str());
	if(!OUTDD.good())
	{
		cerr << "open OUT File error: "<<StatDD<<endl;
		return 0;
	}

	ogzstream OUTRR ((StatRR).c_str());
	if(!OUTRR.good())
	{
		cerr << "open OUT File error: "<<StatRR<<endl;
		return 0;
	}


	ogzstream OUTBB ((StatSite).c_str());
	if(!OUTBB.good())
	{
		cerr << "open OUT File error: "<<StatSite<<endl;
		return 0;
	}

	ogzstream OUTCC ((StatBlocks).c_str());
	if(!OUTCC.good())
	{
		cerr << "open OUT File error: "<<StatBlocks<<endl;
		return 0;
	}
	OUTDD<<"#Region PairWise D'"<<endl;
	OUTRR<<"#Region PairWise R^2"<<endl;

	vector<string> inf;
	split((paraFA04->Region),inf,":-");
	llong StartSite=atoi(inf[1].c_str());
	llong EndSite=atoi(inf[2].c_str());

	string chrName=inf[0];
	unsigned short int ValueE=2;
	map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(chrName);
	if (itSNP == SNPList.end())
	{
		cerr<<"\t\tInPut Para -Region  chromosome ["<<chrName<<"]  can't be found in the SNP dataset\n";		
		return 0;
	}  

	llong BeginRegion=StartSite;   if (BeginRegion<0) { BeginRegion=0;}	
	llong EndRegion=EndSite;

	map<llong, vector <BaseType> >  :: iterator key2=(itSNP->second).begin(); ;

	map <llong, vector <BaseType> >  NewSNPList ;

	int count =0;
	int TotalSNP=0;

	statementVar Var;
	Var.Asize= (key2->second).size();

	for(  ; key2!=(itSNP->second).end(); key2++)
	{
		if ( ( key2->first )<  BeginRegion)
		{
			continue;
		}
		if  ( ( key2->first )> EndRegion)
		{
			break;
		}
		TotalSNP++;
		NewSNPList.insert(map <llong, vector <BaseType> >  :: value_type(key2->first,(key2->second))) ;
	}

	cout<<"##Start Region Cal...\t :"<<chrName<<" "<<BeginRegion<<" "<<EndRegion<<"; In This Region TotalSNP Number is "<<TotalSNP<<endl;
	if (TotalSNP < 3)
	{
		cerr<<"LDBlocks should be with Region SNP Number  > 3 \n";
		return 0;
	}
	else if ( TotalSNP > 168888)
	{
		cerr<<"Warning: LDBlocks Region SNP Number too much,you may use the small region or more stringent conditions to filter the SNP\n";
	}

	//TotalSNP++;
	double **ArrDD = new double *[TotalSNP];
	llong *ArrSite = new llong [TotalSNP];
	bool *TFSite = new bool [TotalSNP];

	for( int i = 0; i < TotalSNP; i++ )
	{
		ArrDD[i] = new double [TotalSNP];
		TFSite[i]=false;
		ArrSite[i]=0;
		for (int j =0 ; j< TotalSNP ; j++)
		{
			ArrDD[i][j]=1.00;
		}
	}



	double ResultRR;
	double ResultDD;
	key2=NewSNPList.begin();
	map<llong, vector <BaseType> >  :: iterator key2_se ;
	int Count=0;

	for(  ; key2!=NewSNPList.end(); key2++ )
	{
		OUTBB<<chrName<<"\t"<<(key2->first)<<endl;
		ArrSite[Count]=(key2->first);

		key2_se=key2  ;  key2_se++ ;
		int SedCount=Count+1;
		if ( key2_se!=NewSNPList.end() )
		{

			cal_RR_D_Both(key2->second , key2_se->second ,ResultRR,ResultDD, Var);

			OUTDD<<setprecision(3)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<ResultDD;
			OUTRR<<setprecision(3)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<ResultRR;

			ArrDD[Count][SedCount]=ResultDD;
			ArrDD[SedCount][Count]=ResultDD;

			key2_se++ ;
			SedCount++;
			for(  ; key2_se!=NewSNPList.end();  key2_se++)
			{
				cal_RR_D_Both(key2->second , key2_se->second ,ResultRR,ResultDD, Var);
				OUTDD<<"\t"<<setprecision(3)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<ResultDD;
				OUTRR<<"\t"<<setprecision(3)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<ResultRR;
				ArrDD[Count][SedCount]=ResultDD;
				ArrDD[SedCount][Count]=ResultDD;
				SedCount++;
			}
			OUTDD<<endl;
			OUTRR<<endl;
		}
		Count++;
	}

	OUTDD.close();
	OUTBB.close();
	OUTRR.close();

	//cerr<<paraFA04->MinRR<<"\t"<<paraFA04->TagRR<<endl;
	OUTCC<<"#chr\tStart\tEnd\tSNPNumber\tTagSNPList\n";
	cerr<<"find blocks...\n";
	if (  (paraFA04->Method)==2 )
	{
		int  SiteStart=0;
		int  SiteEnd=TotalSNP-1;
		while ( SiteStart < TotalSNP )
		{
			for ( SiteEnd=TotalSNP-1   ;  SiteStart<SiteEnd ;SiteEnd--)
			{
				int CountALL=0;
				int CountHLD=0;
				for (int ii=SiteStart+1  ; ii<=SiteEnd  ; ii++)
				{
					CountALL++;
					if (ArrDD[SiteStart][ii]>(paraFA04->MinRR))
					{
						CountHLD++;
					}
				}

				for (int ii=SiteEnd-1  ; ii>SiteStart  ; ii--)
				{
					CountALL++;
					if (ArrDD[SiteStart][SiteEnd]>(paraFA04->MinRR))
					{
						CountHLD++;
					}
				}

				if ( CountHLD >= (CountALL*(paraFA04->RatioRR)))
				{
					break;
				}
			}

			int Length=ArrSite[SiteEnd]-ArrSite[SiteStart];
			if (Length>(paraFA04->SBMinLen) )
			{
				OUTCC<<chrName<<"\t"<<ArrSite[SiteStart]<<"\t"<<ArrSite[SiteEnd]<<"\t"<<SiteEnd-SiteStart+1<<"\t"<<ArrSite[SiteStart];
				TFSite[SiteStart]=true;
				for (int RREE=SiteStart+1; RREE<=SiteEnd; RREE++)
				{
					bool TrunF=true ;
					for (int kk=SiteStart ; kk<RREE ; kk++)
					{
						if ( (TFSite[kk])  &&  (ArrDD[kk][RREE]>(paraFA04->TagRR)) )
						{
							TrunF=false;
						}
					}
					TFSite[RREE]=TrunF;
					if (TrunF){OUTCC<<","<<ArrSite[RREE];}
				}
				OUTCC<<"\n";
			}
			SiteStart=SiteEnd+1;
		}
	}
	else if  ( (paraFA04->Method)==3 )
	{
		int  SiteStart=0;
		int  SiteEnd=TotalSNP-1;
		while ( SiteStart < TotalSNP )
		{
			for ( SiteEnd=TotalSNP-1   ;  SiteStart<SiteEnd ;SiteEnd--)
			{
				int CountALL=0;
				int CountHLD=0;
				for (int ii=SiteStart+1  ; ii<=SiteEnd  ; ii++)
				{
					CountALL++;
					if (ArrDD[ii][ii-1]>(paraFA04->MinRR))
					{
						CountHLD++;
					}
				}
				if ( CountHLD >= (CountALL*(paraFA04->RatioRR)))
				{
					break;
				}
			}

			int Length=ArrSite[SiteEnd]-ArrSite[SiteStart];
			if (Length>(paraFA04->SBMinLen) )
			{
				OUTCC<<chrName<<"\t"<<ArrSite[SiteStart]<<"\t"<<ArrSite[SiteEnd]<<"\t"<<SiteEnd-SiteStart+1<<"\t"<<ArrSite[SiteStart];
				TFSite[SiteStart]=true;
				for (int RREE=SiteStart+1; RREE<=SiteEnd; RREE++)
				{
					bool TrunF=true ;
					for (int kk=SiteStart ; kk<RREE ; kk++)
					{
						if ( (TFSite[kk])  &&  (ArrDD[kk][RREE]>(paraFA04->TagRR)) )
						{
							TrunF=false;
						}
					}
					TFSite[RREE]=TrunF;
					if (TrunF){OUTCC<<","<<ArrSite[RREE];}
				}
				OUTCC<<"\n";
			}
			SiteStart=SiteEnd+1;
			//			cerr<<Length<<"\t"<<ArrDD[ii][ii-1]<<"\t"<<ii<<endl;
		}
	}

	OUTCC.close();

	for(int i = 0; i <TotalSNP; i++)
	{
		delete[] ArrDD[i];
	}
	delete[] ArrDD;
	delete[] ArrSite;
	delete[] TFSite;


	if  ( (paraFA04->Method)==1 )
	{
		stringstream   sstrmMAF , sstrmMiss;
		sstrmMAF  <<  (para_18->MAF)  ;
		sstrmMiss <<  (para_18->Miss) ;
		string Plink=(paraFA04->BinDir)+"/plink";
		if ( access(Plink.c_str(), 0) == 0 )
		{
			Plink=Plink+" --vcf "+(paraFA04->InStr1)+" --geno  " + sstrmMiss.str() +"  -maf " + sstrmMAF.str()  + " --blocks   no-pheno-req  --allow-extra-chr   --out   "+ (paraFA04->InStr2)+".plink > /dev/null  ";
			std::system(Plink.c_str());
			Plink="rm -rf "+(paraFA04->InStr2)+".blocks.gz  ;  mv  "+(paraFA04->InStr2)+".plink.blocks.det  "+(paraFA04->InStr2)+".blocks  ; gzip  " + (paraFA04->InStr2)+".blocks ; rm -rf  "+   (paraFA04->InStr2)+".plink*  " ;
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
				cout <<"\twarning: can't find the [plink] in your $PATH ; no Plink Blocks Figure Out"<<endl;
			}
			else
			{
			  Plink=Plink+" --vcf  "+(paraFA04->InStr1)+" --geno  " + sstrmMiss.str() +"  -maf " + sstrmMAF.str()  + "  --blocks    no-pheno-req  --allow-extra-chr     --out  "+ (paraFA04->InStr2)+".plink > /dev/null  ";
			  std::system(Plink.c_str());
			  Plink="rm -rf "+(paraFA04->InStr2)+".blocks.gz ;  mv  "+(paraFA04->InStr2)+".plink.blocks.det  "+(paraFA04->InStr2)+".blocks  ; gzip  " + (paraFA04->InStr2)+".blocks ;  rm  -rf " +   (paraFA04->InStr2)+".plink*  "  ;
			  std::system(Plink.c_str());
			}
		}
	}




	return 1;

}










int Region_PairWise_BothB(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > >  & SNPList, int & Flag_for_pro)
{
	string StatDD  =(paraFA04->InStr2)+".TriangleB.gz";
	string StatRR  =(paraFA04->InStr2)+".TriangleV.gz";
	string StatSite=(paraFA04->InStr2)+".site.gz";
	string StatBlocks=(paraFA04->InStr2)+".blocks.gz";
	ogzstream OUTRR ((StatRR).c_str());
	if(!OUTRR.good())
	{
		cerr << "open OUT File error: "<<StatRR<<endl;
		return 0;
	}

	ogzstream OUTDD ((StatDD).c_str());
	if(!OUTDD.good())
	{
		cerr << "open OUT File error: "<<StatDD<<endl;
		return 0;
	}



	ogzstream OUTBB ((StatSite).c_str());
	if(!OUTBB.good())
	{
		cerr << "open OUT File error: "<<StatSite<<endl;
		return 0;
	}

	ogzstream OUTCC ((StatBlocks).c_str());
	if(!OUTCC.good())
	{
		cerr << "open OUT File error: "<<StatBlocks<<endl;
		return 0;
	}
	OUTRR<<"#Region PairWise R^2"<<endl;
	OUTDD<<"#Region PairWise D'"<<endl;

	vector<string> inf;
	split((paraFA04->Region),inf,":-");
	llong StartSite=atoi(inf[1].c_str());
	llong EndSite=atoi(inf[2].c_str());

	string chrName=inf[0];
	unsigned short int ValueE=2;
	map <string,map <llong, vector <BaseType> > >  :: iterator itSNP=SNPList.find(chrName);
	if (itSNP == SNPList.end())
	{
		cerr<<"\t\tInPut Para -Region  chromosome ["<<chrName<<"]  can't be found in the SNP dataset\n";		
		return 0;
	}  

	llong BeginRegion=StartSite;   if (BeginRegion<0) { BeginRegion=0;}	
	llong EndRegion=EndSite;

	map<llong, vector <BaseType> >  :: iterator key2=(itSNP->second).begin(); ;

	map <llong, vector <BaseType> >  NewSNPList ;

	int count =0;
	int TotalSNP=0;

	statementVar Var;
	Var.Asize= (key2->second).size();

	for(  ; key2!=(itSNP->second).end(); key2++)
	{
		if ( ( key2->first )<  BeginRegion)
		{
			continue;
		}
		if  ( ( key2->first )> EndRegion)
		{
			break;
		}
		TotalSNP++;
		NewSNPList.insert(map <llong, vector <BaseType> >  :: value_type(key2->first,(key2->second))) ;
	}

	cout<<"##Start Region Cal...\t :"<<chrName<<" "<<BeginRegion<<" "<<EndRegion<<"; In This Region TotalSNP Number is "<<TotalSNP<<endl;
	if (TotalSNP < 3)
	{
		cerr<<"LDBlocks should be with Region SNP Number  > 3 \n";
		return 0;
	}
	else if ( TotalSNP > 168888)
	{
		cerr<<"Warning: LDBlocks Region SNP Number too much,you may use the small region or more stringent conditions to filter the SNP\n";
	}

	//TotalSNP++;
	double **ArrDD = new double *[TotalSNP];
	llong *ArrSite = new llong [TotalSNP];
	bool *TFSite = new bool [TotalSNP];

	for( int i = 0; i < TotalSNP; i++ )
	{
		ArrDD[i] = new double [TotalSNP];
		TFSite[i]=false;
		ArrSite[i]=0;
		for (int j =0 ; j< TotalSNP ; j++)
		{
			ArrDD[i][j]=1.00;
		}
	}



	double ResultRR;
	double ResultDD;
	key2=NewSNPList.begin();
	map<llong, vector <BaseType> >  :: iterator key2_se ;
	int Count=0;

	for(  ; key2!=NewSNPList.end(); key2++ )
	{
		OUTBB<<chrName<<"\t"<<(key2->first)<<endl;
		ArrSite[Count]=(key2->first);

		key2_se=key2  ;  key2_se++ ;
		int SedCount=Count+1;
		if ( key2_se!=NewSNPList.end() )
		{

			cal_RR_D_Both(key2->second , key2_se->second ,ResultRR,ResultDD, Var);

			OUTRR<<setprecision(3)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<ResultRR;
			OUTDD<<setprecision(3)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<ResultDD;

			ArrDD[Count][SedCount]=ResultRR;
			ArrDD[SedCount][Count]=ResultRR;

			key2_se++ ;
			SedCount++;
			for(  ; key2_se!=NewSNPList.end();  key2_se++)
			{
				cal_RR_D_Both(key2->second , key2_se->second ,ResultRR,ResultDD, Var);
				OUTRR<<"\t"<<setprecision(3)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<ResultRR;
				OUTDD<<"\t"<<setprecision(3)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<ResultDD;
				ArrDD[Count][SedCount]=ResultRR;
				ArrDD[SedCount][Count]=ResultRR;
				SedCount++;
			}
			OUTRR<<endl;
			OUTDD<<endl;
		}
		Count++;
	}

	OUTBB.close();
	OUTRR.close();
	OUTDD.close();

	//cerr<<paraFA04->MinRR<<"\t"<<paraFA04->TagRR<<endl;
	OUTCC<<"#chr\tStart\tEnd\tSNPNumber\tTagSNPList\n";
	cerr<<"find blocks...\n";
	if (  (paraFA04->Method)==2 )
	{
		int  SiteStart=0;
		int  SiteEnd=TotalSNP-1;
		while ( SiteStart < TotalSNP )
		{
			for ( SiteEnd=TotalSNP-1   ;  SiteStart<SiteEnd ;SiteEnd--)
			{
				int CountALL=0;
				int CountHLD=0;
				for (int ii=SiteStart+1  ; ii<=SiteEnd  ; ii++)
				{
					CountALL++;
					if (ArrDD[SiteStart][ii]>(paraFA04->MinRR))
					{
						CountHLD++;
					}
				}

				for (int ii=SiteEnd-1  ; ii>SiteStart  ; ii--)
				{
					CountALL++;
					if (ArrDD[SiteStart][SiteEnd]>(paraFA04->MinRR))
					{
						CountHLD++;
					}
				}

				if ( CountHLD >= (CountALL*(paraFA04->RatioRR)))
				{
					break;
				}
			}

			int Length=ArrSite[SiteEnd]-ArrSite[SiteStart];
			if (Length>(paraFA04->SBMinLen) )
			{
				OUTCC<<chrName<<"\t"<<ArrSite[SiteStart]<<"\t"<<ArrSite[SiteEnd]<<"\t"<<SiteEnd-SiteStart+1<<"\t"<<ArrSite[SiteStart];
				TFSite[SiteStart]=true;
				for (int RREE=SiteStart+1; RREE<=SiteEnd; RREE++)
				{
					bool TrunF=true ;
					for (int kk=SiteStart ; kk<RREE ; kk++)
					{
						if ( (TFSite[kk])  &&  (ArrDD[kk][RREE]>(paraFA04->TagRR)) )
						{
							TrunF=false;
						}
					}
					TFSite[RREE]=TrunF;
					if (TrunF){OUTCC<<","<<ArrSite[RREE];}
				}
				OUTCC<<"\n";
			}
			SiteStart=SiteEnd+1;
		}
	}
	else if  ( (paraFA04->Method)==3 )
	{
		int  SiteStart=0;
		int  SiteEnd=TotalSNP-1;
		while ( SiteStart < TotalSNP )
		{
			for ( SiteEnd=TotalSNP-1   ;  SiteStart<SiteEnd ;SiteEnd--)
			{
				int CountALL=0;
				int CountHLD=0;
				for (int ii=SiteStart+1  ; ii<=SiteEnd  ; ii++)
				{
					CountALL++;
					if (ArrDD[ii][ii-1]>(paraFA04->MinRR))
					{
						CountHLD++;
					}
				}
				if ( CountHLD >= (CountALL*(paraFA04->RatioRR)))
				{
					break;
				}
			}

			int Length=ArrSite[SiteEnd]-ArrSite[SiteStart];
			if (Length>(paraFA04->SBMinLen) )
			{
				OUTCC<<chrName<<"\t"<<ArrSite[SiteStart]<<"\t"<<ArrSite[SiteEnd]<<"\t"<<SiteEnd-SiteStart+1<<"\t"<<ArrSite[SiteStart];
				TFSite[SiteStart]=true;
				for (int RREE=SiteStart+1; RREE<=SiteEnd; RREE++)
				{
					bool TrunF=true ;
					for (int kk=SiteStart ; kk<RREE ; kk++)
					{
						if ( (TFSite[kk])  &&  (ArrDD[kk][RREE]>(paraFA04->TagRR)) )
						{
							TrunF=false;
						}
					}
					TFSite[RREE]=TrunF;
					if (TrunF){OUTCC<<","<<ArrSite[RREE];}
				}
				OUTCC<<"\n";
			}
			SiteStart=SiteEnd+1;
			//			cerr<<Length<<"\t"<<ArrDD[ii][ii-1]<<"\t"<<ii<<endl;
		}
	}

	OUTCC.close();

	for(int i = 0; i <TotalSNP; i++)
	{
		delete[] ArrDD[i];
	}
	delete[] ArrDD;
	delete[] ArrSite;
	delete[] TFSite;


	if  ( (paraFA04->Method)==1 )
	{
		stringstream   sstrmMAF , sstrmMiss;
		sstrmMAF  <<  (para_18->MAF)  ;
		sstrmMiss <<  (para_18->Miss) ;
		string Plink=(paraFA04->BinDir)+"/plink";
		if ( access(Plink.c_str(), 0) == 0 )
		{
			Plink=Plink+" --vcf "+(paraFA04->InStr1)+" --geno  " + sstrmMiss.str() +"  -maf " + sstrmMAF.str()  + " --blocks   no-pheno-req  --allow-extra-chr   --out   "+ (paraFA04->InStr2)+".plink > /dev/null  ";
			std::system(Plink.c_str());
			Plink="rm -rf "+(paraFA04->InStr2)+".blocks.gz  ;  mv  "+(paraFA04->InStr2)+".plink.blocks.det  "+(paraFA04->InStr2)+".blocks  ; gzip  " + (paraFA04->InStr2)+".blocks ; rm -rf  "+   (paraFA04->InStr2)+".plink*  " ;
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
				cout <<"\twarning: can't find the [plink] in your $PATH ; no Plink Blocks Figure Out"<<endl;
			}
			else
			{
			  Plink=Plink+" --vcf  "+(paraFA04->InStr1)+" --geno  " + sstrmMiss.str() +"  -maf " + sstrmMAF.str()  + "  --blocks    no-pheno-req  --allow-extra-chr     --out  "+ (paraFA04->InStr2)+".plink > /dev/null  ";
			  std::system(Plink.c_str());
			  Plink="rm -rf "+(paraFA04->InStr2)+".blocks.gz ;  mv  "+(paraFA04->InStr2)+".plink.blocks.det  "+(paraFA04->InStr2)+".blocks  ; gzip  " + (paraFA04->InStr2)+".blocks ;  rm  -rf " +   (paraFA04->InStr2)+".plink*  "  ;
			  std::system(Plink.c_str());
			}
		}
	}




	return 1;

}










#endif
