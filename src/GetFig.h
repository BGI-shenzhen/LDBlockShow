#ifndef GetFig_H_
#define GetFig_H_


using namespace std;

int GetFig_Fun(In3str1v * paraFA04, Para_18 * para_18 )
{
	string binPath=(paraFA04->BinDir)+"/ShowLDSVG";
	string command=binPath+"\t-InPreFix\t"+(paraFA04->InStr2)+"\t-OutPut\t"+(paraFA04->InStr2)+".svg";
	if  (!(paraFA04->InGWAS).empty())
	{
		command=command+"\t"+(paraFA04->InGWAS);
	}

	if ( access(binPath.c_str(), 0) == 0 )
	{		
		std::system(command.c_str()) ;
	}	
	else
	{
		char   buf[2048]={'\0'};
		string cc="which  ShowLDSVG  2> /dev/null ";
		memset( buf, '\0', sizeof(buf) );
		FILE   *stream ;
		stream=popen(cc.c_str(),"r") ;
		fread( buf, sizeof(char), sizeof(buf), stream);
		binPath=buf;
		binPath=binPath.substr(0,binPath.length()-1);

		if (binPath == "")
		{
			cout <<"\twarning: can't find the [ShowLDSVG] in your $PATH ; no png Figure Out"<<endl;
			cout <<"\t\tShowLDSVG\t"<<command<<endl;
		}
		else
		{
			std::system(command.c_str()) ;
		}
	}

	return 1;
}





#endif

