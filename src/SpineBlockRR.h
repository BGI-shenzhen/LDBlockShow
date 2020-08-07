#ifndef SpineBlockRR
#define SpineBlockRR

using namespace std;


int SpineBlockRRFun(In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > >  & SNPList, int & Flag_for_pro)
{

	string StatBlocks=(paraFA04->InStr2)+".blocks.gz";	
	ogzstream OUTCC ((StatBlocks).c_str());
	if(!OUTCC.good())
	{
		cerr << "open OUT File error: "<<StatBlocks<<endl;
		return 0;
	}

	OUTCC<<"#chr\tStart\tEnd\tSNPNumber\tTagSNPList\n";


	map <string,map <llong, vector <BaseType> > >  :: iterator itSNPChr;
	map <llong, vector <BaseType> >  :: iterator  key2 ;
	map <llong, vector <BaseType> >  :: iterator key2_se ;

	double CalResult;
	statementVar Var;
	map  <llong,bool > TagSNPSite;

	for(  itSNPChr=SNPList.begin() ; itSNPChr!=SNPList.end(); itSNPChr++)
	{

		key2=(itSNPChr->second).begin();
		Var.Asize= (key2->second).size();
		int  Length=0;
		llong  SiteStart=0;
		llong  SiteEnd=0;


		for(  ; key2!=(itSNPChr->second).end(); key2++ )
		{
			key2_se=key2  ;  key2_se++ ;

			if ( key2_se!=(itSNPChr->second).end() )
			{
				cal_RR_MA(key2->second , key2_se->second ,CalResult, Var);
				if (CalResult>(paraFA04->MinRR))
				{
					if  (Length==0)	
					{
						SiteStart=key2->first;
						SiteEnd=key2_se->first;
						Length=1;
						TagSNPSite.clear();
						TagSNPSite.insert(map <llong, bool >  :: value_type(SiteStart,true)) ;
						if  ( CalResult<(paraFA04->TagRR) )
						{
							TagSNPSite.insert(map <llong, bool >  :: value_type(SiteEnd,true)) ;
						}
					}
					else
					{
						SiteEnd=key2_se->first;
						if  ( CalResult<(paraFA04->TagRR) )
						{
							TagSNPSite.insert(map <llong, bool >  :: value_type(SiteEnd,true)) ;
						}
						Length++;
					}
				}
				else
				{
					if  (Length>2)
					{
						OUTCC<<itSNPChr->first<<"\t"<<SiteStart<<"\t"<<SiteEnd<<"\t"<<Length+1;
						map<llong, bool>::iterator TagSNPit= TagSNPSite.begin();
						OUTCC<<"\t"<<TagSNPit->first;
						TagSNPit++;
						for (; TagSNPit != TagSNPSite.end(); TagSNPit++)
						{
							OUTCC<<","<<TagSNPit->first;
						}			
						OUTCC<<"\n";
						Length=0;
					}
				}
			}
		}
	}
		OUTCC.close();
		return 1;

}




#endif
