#ifndef comm_H_
#define comm_H_


using namespace std;
typedef int  llong ;


//////////////////////////////// some small funtion /////////////////////////////////////////////

bool GetBinPath ( std::string path  ,std::string& BinPath )
{
	char   buf[2048]={'\0'};
	BinPath.clear() ;
	if ( path.rfind('/') == string::npos)
	{
		string cc="which  "+path+"   2> /dev/null ";
		FILE   *stream ;
		//cout <<cc<<endl;
		stream=popen(cc.c_str(),"r") ;
		fread( buf, sizeof(char), sizeof(buf), stream);
		BinPath=buf;
		//cout <<"aa\t"<<BinPath<<"\teee"<<endl;
		BinPath=BinPath.substr(0,BinPath.rfind('/'));
		//cout <<"aa\t"<<BinPath<<"\teee"<<endl;
	}
	else
	{
		if (path[0]  == '/' )
		{
			BinPath=path.substr(0, path.rfind('/')+1);
		}
		else if   (path[0]  == '.' )
		{

			if(!getcwd(buf, 1000))
			{
				return false ;
			}
			BinPath = buf ;

			if (path[1]  == '/')
			{
				string cutOFF =path.substr(1, path.rfind('/'));
				BinPath = BinPath+cutOFF;
			}
			else
			{
				string cutOFF =path.substr(0, path.rfind('/'));
				BinPath = BinPath+"/"+cutOFF;
			}
		}
		else
		{
			string cutOFF =path.substr(0, path.rfind('/')+1);
			if(!getcwd(buf, 1000))
			{
				return false ;
			}
			BinPath = buf ;
			BinPath = BinPath+"/"+cutOFF;
		}
	}
	return true ;
}

string Int2Str (int A )
{
	stringstream   sstrm ;
	sstrm  <<  A ;
	return  sstrm.str();
}

/*
string int2str( unsigned short int val )
{
	ostringstream out;
	out<<val;
	return out.str();
}
*/

inline void  LogLackArg( string  flag )
{
	cerr << "\t\tLack Argument for [ -"<<flag<<" ]"<<endl;
}

/*
   inline string add_Asuffix ( string path )
   {
   string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);
   if (ext != "gz")
   {
   path=path+".gz" ;
   }
   return path ;
   }
   */

inline string &  replace_all(string &  str,const  string &  old_Avalue,const string &  new_Avalue)
{
	while(true)   {
		string::size_type  pos(0);
		if(   (pos=str.find(old_Avalue))!=string::npos   )
			str.replace(pos,old_Avalue.length(),new_Avalue);
		else   break;
	}
	return   str;
}


inline void split(const string& str,vector<string>& tokens,  const string& delimiters = " ")
{
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	string::size_type pos     = str.find_first_of(delimiters, lastPos);
	while (string::npos != pos || string::npos != lastPos)
	{
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		lastPos = str.find_first_not_of(delimiters, pos);
		pos = str.find_first_of(delimiters, lastPos);
	}
}



bool  GzipTrue( string filepath )
{
	FILE   *infile;
	unsigned char buf[3] ={'\0'};
	infile = fopen(filepath.c_str(), "rb");
	int rc = fread(buf,sizeof(unsigned char), 2,infile);
	fclose(infile);
	if  (buf[0] == 31 &&  buf[1] == 139)
	{
		return true ;
	}
	else
	{
		return false ;
	}
}



///////////////////swimming in the sky & flying in the sea/////////////////////////////





#endif // comm_H_  ;


///////////////////



////////////////////////swimming in the sea & flying in the sky //////////////////


