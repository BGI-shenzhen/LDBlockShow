#ifndef calculate_H_
#define calculate_H_


using namespace std;


/////////////////////////       Cal method1  RR  /////////////////////////////////
int cal_RR_MA ( vector<BaseType>  & Base1  , vector<BaseType> &  Base2 ,  double &  CalResult, statementVar & Var )
{
	Var.DDE[0][0]=0;	Var.DDE[0][1]=0;	Var.DDE[1][0]=0;	Var.DDE[1][1]=0;
	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Var.DDE[(Base1[Var.i].Value)][(Base2[Var.i].Value)]++ ;
	}
	
	Var.tmpAA=(Var.DDE[1][1])+(Var.DDE[1][0]);
	if (Var.tmpAA==0)
	{
		return 0 ;
	}
	if ( (Var.DDE[1][1]+Var.DDE[0][1])==0)
	{
		return 0 ;
	}

	Var.ALL_count=Var.DDE[0][0]+Var.DDE[0][1]+Var.tmpAA;
	Var.probHaps[0]=((Var.DDE[0][0])/Var.ALL_count);
	Var.probHaps[1]=((Var.DDE[0][1])/Var.ALL_count);
	Var.probHaps[2]=((Var.DDE[1][0])/Var.ALL_count);

	Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
	Var.pA2 = Var.probHaps[0]+Var.probHaps[2];

	Var.Cal_B=(Var.pA1)*(Var.pA2);
	Var.Cal_A  = 1.0-(Var.pA1+Var.pA2)+Var.Cal_B ;

	if  (Var.Cal_A==0  || Var.Cal_B==0 )
	{
		if (Var.probHaps[0] < 1e-10) { Var.probHaps[0]=1e-10;}
		if (Var.probHaps[1] < 1e-10) { Var.probHaps[1]=1e-10;}
		if (Var.probHaps[2] < 1e-10) { Var.probHaps[2]=1e-10;}

		Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
		Var.pA2 = Var.probHaps[0]+Var.probHaps[2];

		Var.Cal_B = (Var.pA1)*(Var.pA2);
	    Var.Cal_A  = 1.0-(Var.pA1+Var.pA2)+Var.Cal_B ;

	}

	Var.D_A = Var.probHaps[0]-Var.Cal_B ;
	CalResult = (Var.D_A*Var.D_A)/(Var.Cal_A*Var.Cal_B);

	return 1;
}





/////////////////////////       Cal method1  D'   /////////////////////////////////

 int cal_RR_D(vector<BaseType>  & Base1  , vector<BaseType> &  Base2 ,  double  &  CalResult, statementVar & Var)
{

	Var.DDE[0][0]=0;	Var.DDE[0][1]=0;	Var.DDE[1][0]=0;	Var.DDE[1][1]=0;

	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Var.DDE[(Base1[Var.i].Value)][(Base2[Var.i].Value)]++;
	}

	Var.tmpAA=Var.DDE[1][1]+Var.DDE[1][0];

	if (Var.tmpAA==0)
	{
		return 0 ;
	}
	if ( (Var.DDE[1][1]+Var.DDE[0][1])==0)
	{
		return 0 ;
	}

	Var.ALL_count=Var.DDE[0][0]+Var.DDE[0][1]+Var.tmpAA;
	Var.probHaps[0]=((Var.DDE[0][0])/Var.ALL_count);
	Var.probHaps[1]=((Var.DDE[0][1])/Var.ALL_count);
	Var.probHaps[2]=((Var.DDE[1][0])/Var.ALL_count);

	Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
	Var.pB1 = 1.0-Var.pA1;
	Var.pA2 = Var.probHaps[0]+Var.probHaps[2];
	Var.pB2 = 1.0-Var.pA2;

	Var.XpA1_pA2=Var.pA1*Var.pA2;
	Var.D_A = Var.probHaps[0]-Var.XpA1_pA2 ;

	if  (Var.D_A>0)
	{
		Var.Cal_A=Var.pB1*Var.pA2;
		Var.Cal_B=Var.pA1*Var.pB2; 
	}
	else 
	{
		Var.D_A = 0.0-Var.D_A;
		Var.Cal_A = (Var.pB1)*(Var.pB2);
		Var.Cal_B = Var.XpA1_pA2 ;
	}

	Var.D_max=Var.Cal_A;	if  (Var.Cal_A>Var.Cal_B)	{		Var.D_max=Var.Cal_B;	}

	if  (Var.D_max==0)
	{
		if (Var.probHaps[0] < 1e-10) { Var.probHaps[0]=1e-10;}
		if (Var.probHaps[1] < 1e-10) { Var.probHaps[1]=1e-10;}
		if (Var.probHaps[2] < 1e-10) { Var.probHaps[2]=1e-10;}

		Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
		Var.pB1 = 1.0-Var.pA1;
		Var.pA2 = Var.probHaps[0]+Var.probHaps[2];
		Var.pB2 = 1.0-Var.pA2;

		Var.XpA1_pA2=Var.pA1*Var.pA2;
		Var.D_A = Var.probHaps[0]-Var.XpA1_pA2 ;

		if  (Var.D_A>0)
		{
			Var.Cal_A=Var.pB1*Var.pA2;
			Var.Cal_B=Var.pA1*Var.pB2; 
		}
		else 
		{
			Var.D_A = 0.0-Var.D_A;
			Var.Cal_A = (Var.pB1)*(Var.pB2);
			Var.Cal_B = Var.XpA1_pA2 ;
		}
		Var.D_max=Var.Cal_A;	if  (Var.Cal_A>Var.Cal_B)	{Var.D_max=Var.Cal_B;}
	}

	CalResult  = Var.D_A/Var.D_max;
//	CalResult.D = Var.D_A/Var.D_max;
//	CalResult.RR = (Var.D_A/Var.Cal_A)*(Var.D_A/Var.Cal_B);

	return 1;
}





 int cal_RR_D_Both(vector<BaseType>  & Base1  , vector<BaseType> &  Base2 ,  double  &  Result_RR, double  &  Result_DD ,  statementVar & Var)
{

	Var.DDE[0][0]=0;	Var.DDE[0][1]=0;	Var.DDE[1][0]=0;	Var.DDE[1][1]=0;

	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Var.DDE[(Base1[Var.i].Value)][(Base2[Var.i].Value)]++;
	}

	Var.tmpAA=Var.DDE[1][1]+Var.DDE[1][0];

	if (Var.tmpAA==0)
	{
		return 0 ;
	}
	if ( (Var.DDE[1][1]+Var.DDE[0][1])==0)
	{
		return 0 ;
	}

	Var.ALL_count=Var.DDE[0][0]+Var.DDE[0][1]+Var.tmpAA;
	Var.probHaps[0]=((Var.DDE[0][0])/Var.ALL_count);
	Var.probHaps[1]=((Var.DDE[0][1])/Var.ALL_count);
	Var.probHaps[2]=((Var.DDE[1][0])/Var.ALL_count);

	Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
	Var.pB1 = 1.0-Var.pA1;
	Var.pA2 = Var.probHaps[0]+Var.probHaps[2];
	Var.pB2 = 1.0-Var.pA2;

	Var.XpA1_pA2=Var.pA1*Var.pA2;
	Var.D_A = Var.probHaps[0]-Var.XpA1_pA2 ;

	if  (Var.D_A>0)
	{
		Var.Cal_A=Var.pB1*Var.pA2;
		Var.Cal_B=Var.pA1*Var.pB2; 
	}
	else 
	{
		Var.D_A = 0.0-Var.D_A;
		Var.Cal_A = (Var.pB1)*(Var.pB2);
		Var.Cal_B = Var.XpA1_pA2 ;
	}

	Var.D_max=Var.Cal_A;	if  (Var.Cal_A>Var.Cal_B)	{		Var.D_max=Var.Cal_B;	}

	if  (Var.D_max==0)
	{
		if (Var.probHaps[0] < 1e-10) { Var.probHaps[0]=1e-10;}
		if (Var.probHaps[1] < 1e-10) { Var.probHaps[1]=1e-10;}
		if (Var.probHaps[2] < 1e-10) { Var.probHaps[2]=1e-10;}

		Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
		Var.pB1 = 1.0-Var.pA1;
		Var.pA2 = Var.probHaps[0]+Var.probHaps[2];
		Var.pB2 = 1.0-Var.pA2;

		Var.XpA1_pA2=Var.pA1*Var.pA2;
		Var.D_A = Var.probHaps[0]-Var.XpA1_pA2 ;

		if  (Var.D_A>0)
		{
			Var.Cal_A=Var.pB1*Var.pA2;
			Var.Cal_B=Var.pA1*Var.pB2; 
		}
		else 
		{
			Var.D_A = 0.0-Var.D_A;
			Var.Cal_A = (Var.pB1)*(Var.pB2);
			Var.Cal_B = Var.XpA1_pA2 ;
		}
		Var.D_max=Var.Cal_A;	if  (Var.Cal_A>Var.Cal_B)	{Var.D_max=Var.Cal_B;}
	}

	Result_DD  = Var.D_A/Var.D_max;
	Result_RR = (Var.D_A/Var.Cal_A)*(Var.D_A/Var.Cal_B);

	return 1;
}










#endif // calculate_H_  ;
/////////////////////////////////////////////////////////////////////////////////////////////////

