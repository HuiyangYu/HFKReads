/* The MIT License

   Copyright (c) 2023- by Huiyang Yu, Weiming He, Chunmei Shi.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   */

#ifndef FQ_Comm_H_
#define FQ_Comm_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <thread>
#include <algorithm>
#include <malloc.h>
#include <unordered_map>
#include <map>
#include <cstdio>
#include <zlib.h>
#include "./gzstream.c"
#include "./kseq.h"
#include "./DataClass.h"

//

using namespace std;


KSEQ_INIT(gzFile, gzread)





string &  replace_all(string &  str,const  string &  old_Avalue,const string &  new_Avalue)
{
	while(true)   {
		string::size_type  pos(0);
		if(   (pos=str.find(old_Avalue))!=string::npos   )
			str.replace(pos,old_Avalue.length(),new_Avalue);
		else   break;
	}
	return   str;
}



inline void splitID(string& str,  const string& delimiters = " ")
{
	string::size_type pos = str.find_first_of(delimiters,1);
	if (string::npos!=pos)
	{
		str='>'+str.substr(0, pos);
	}
}


inline void  LogLackArg( string  flag )
{
	cerr << "Error: Lack Argument for [ -"<<flag<<" ]"<<endl;
}

int GetShiftQ ( string FQPath, Para_A24 * P2In )
{
	igzstream INFQ ((FQPath).c_str(),ifstream::in);
	if (INFQ.fail())
	{
		cerr << "Error: open File error: "<<(FQPath)<<endl;
		return 64 ;
	}
	string tmp ,tmp2,tmp3,tmp4 ;
	int minQ=50000;
	int maxQ=0;
	getline(INFQ,tmp);

	string::size_type pos = tmp.find_first_of(" \t",1);
	if (string::npos!=pos)
	{
		tmp=tmp.substr(0, pos);
	}

	tmp2=tmp.substr(tmp.size()-2);
	if  (tmp2 =="/1"  || tmp2 =="/2" )
	{
		P2In->IDAdd=true;
	}

	getline(INFQ,tmp2);
	getline(INFQ,tmp3);
	getline(INFQ,tmp4);
	if ((tmp[0] == '>'  )  &&  (tmp3[0] == '>'  ))
	{
		P2In->ReadLength=tmp2.length();
		P2In->N_Number=int((P2In->ReadLength)*(P2In->N_Ration));
		return 0;
	}

	for (int A=1 ; A<1688 && (!INFQ.eof())  ; A++ )
	{
		getline(INFQ,tmp);
		if (tmp.length()<=0)  { continue  ; }

		if(A%4!=0)
		{
			continue;
		}
		string::size_type SeqQLength =tmp.size();
		if ((P2In->ReadLength)<SeqQLength)
		{
			P2In->ReadLength=SeqQLength;
			P2In->N_Number=int((P2In->ReadLength)*(P2In->N_Ration));
		}
		for(int i=0 ; i<SeqQLength ; i++)
		{
			if(minQ>tmp[i])
			{
				minQ=tmp[i];
			}
			if(maxQ<tmp[i])
			{
				maxQ=tmp[i];
			}
		}
	}
	INFQ.close();

	if(minQ >= 33 &&  minQ <= 78  &&  maxQ >= 33 && maxQ <=78 )
	{
		return 33;
	}
	else if (minQ >= 64  &&  minQ <= 108  &&  maxQ >= 64 && maxQ <= 108)
	{
		return 64;
	}
	else
	{
		return 64 ;
	}
}









inline bool FilterReadLowQ( string & Seq , string & Quli ,Para_A24 * P2In)
{
	short int sumQ=0;
	bool Low=true;
	int NN=0;
	string::size_type EEEnd=Seq.size();
	if (  EEEnd< (P2In->HarfReadLength) ) {return false ;}
	for(string::size_type ix=0 ; ix<EEEnd ; ix++)
	{
		NN+=NArry[Seq[ix]] ;	
		Low=Low & LowArry[Quli[ix]];
		sumQ+=Quli[ix];
	}

	if (NN>(P2In->N_Number)) {return false ;}
	if (!Low)
	{
		return false ;
	}

	if  (  (sumQ/EEEnd) < (P2In->AverQ)  )
	{
		return false  ;
	}
	return true ;
}


inline bool FilterReadNNNQ( string & Seq  ,Para_A24 * P2In)
{
	int NN=0;
	string::size_type EEEnd=Seq.size();
	if (  EEEnd< (P2In->HarfReadLength) ) {return false ;}
	for(string::size_type ix=0 ; ix<EEEnd ; ix++)
	{
		NN+=NArry[Seq[ix]] ;	
	}
	if (NN>(P2In->N_Number)) {return false ;}
	return true ;
}

void FilterFQPE(Para_A24 * P2In,  bool  * PASS, int & Start,int &   End ,vector <string> & AAASSS,vector <string> & AAAQQQ,vector <string> & BBBSSS,vector <string> & BBBQQQ)
{
	bool  QPASS_AA;
	bool  QPASS_BB;
	unordered_map <string, bool > :: iterator  MapIt;
	string  Cat ;
	unordered_map <string, bool >   localPCR ;

	for (int ii=Start; ii<End; ii++)
	{
		Cat= AAASSS[ii]+BBBSSS[ii];
		MapIt=localPCR.find(Cat);
		if ( MapIt != localPCR.end() )
		{
			PASS[ii]=false;
			continue ;
		}
		else
		{
			localPCR[Cat]=true;
		}
		QPASS_AA=FilterReadLowQ( AAASSS[ii],AAAQQQ[ii], P2In);
		QPASS_BB=FilterReadLowQ( BBBSSS[ii],BBBQQQ[ii], P2In);
		PASS[ii]=QPASS_AA & QPASS_BB ;
	}
}



void M1_FilterFQPE(Para_A24 * P2In,  bool  * PASSAAA,bool  * PASSBBB, int & Start,int &   End ,vector <string> & AAASSS,vector <string> & AAAQQQ,vector <string> & BBBSSS,vector <string> & BBBQQQ)
{
	unordered_map <string, bool > :: iterator  MapIt;
	string  Cat;
	unordered_map <string, bool >   localPCR ;
	for (int ii=Start; ii<End; ii++)
	{
		Cat= AAASSS[ii]+BBBSSS[ii];
		MapIt=localPCR.find(Cat);
		if ( MapIt != localPCR.end())
		{
			PASSAAA[ii]=false;
			PASSBBB[ii]=false;
			continue ;
		}
		else
		{
			localPCR[Cat]=true;
		}
		PASSAAA[ii]=FilterReadLowQ(AAASSS[ii],AAAQQQ[ii],P2In);
		PASSBBB[ii]=FilterReadLowQ(BBBSSS[ii],BBBQQQ[ii],P2In);
	}
}



void FilterFAPE(Para_A24 * P2In,  bool  * PASS, int & Start,int &   End ,vector <string> & AAASSS,vector <string> & BBBSSS)
{
	bool  QPASS_AA;
	bool  QPASS_BB;
	unordered_map <string, bool > :: iterator  MapIt;
	string  Cat ;
	unordered_map <string, bool >   localPCR ;

	for (int ii=Start; ii<End; ii++)
	{
		Cat= AAASSS[ii]+BBBSSS[ii];
		MapIt=localPCR.find(Cat);
		if ( MapIt != localPCR.end() )
		{
			PASS[ii]=false;
			continue ;
		}
		else
		{
			localPCR[Cat]=true;
		}
		QPASS_AA=FilterReadNNNQ( AAASSS[ii], P2In);
		QPASS_BB=FilterReadNNNQ( BBBSSS[ii], P2In);
		PASS[ii]=QPASS_AA & QPASS_BB ;
	}
}



void M1_FilterFAPE(Para_A24 * P2In,  bool  * PASSAAA, bool  * PASSBBB , int & Start,int &   End ,vector <string> & AAASSS,vector <string> & BBBSSS)
{
	unordered_map <string, bool > :: iterator  MapIt;
	string  Cat ;
	unordered_map <string, bool >   localPCR ;

	for (int ii=Start; ii<End; ii++)
	{
		Cat= AAASSS[ii]+BBBSSS[ii];
		MapIt=localPCR.find(Cat);
		if ( MapIt != localPCR.end() )
		{
			PASSAAA[ii]=false;
			PASSBBB[ii]=false;
			continue ;
		}
		else
		{
			localPCR[Cat]=true;
		}
		PASSAAA[ii]=FilterReadNNNQ( AAASSS[ii], P2In);
		PASSBBB[ii]=FilterReadNNNQ( BBBSSS[ii], P2In);
	}
}



void FilterFQSE(Para_A24 * P2In,  bool  * PASS, int & Start,int &   End ,vector <string> & AAASSS,vector <string> & AAAQQQ)
{
	unordered_map <string, bool > :: iterator  MapIt;
	string  Cat ;
	unordered_map <string, bool >   localPCR ;

	for (int ii=Start; ii<End; ii++)
	{
		Cat= AAASSS[ii];
		MapIt=localPCR.find(Cat);
		if ( MapIt != localPCR.end() )
		{
			PASS[ii]=false;
			continue ;
		}
		else
		{
			localPCR[Cat]=true;
		}
		PASS[ii]=FilterReadLowQ( AAASSS[ii],AAAQQQ[ii], P2In);
	}
}




void FilterFASE(Para_A24 * P2In,  bool  * PASS, int & Start,int &   End ,vector <string> & AAASSS)
{
	unordered_map <string, bool > :: iterator  MapIt;
	string  Cat ;
	unordered_map <string, bool >   localPCR ;

	for (int ii=Start; ii<End; ii++)
	{
		Cat= AAASSS[ii];
		MapIt=localPCR.find(Cat);
		if ( MapIt != localPCR.end() )
		{
			PASS[ii]=false;
			continue ;
		}
		else
		{
			localPCR[Cat]=true;
		}
		PASS[ii]=FilterReadNNNQ( AAASSS[ii], P2In);
	}
}


/*
inline void splitIDM1(string& str,  const string& delimiters = " ")
{
	string::size_type pos = str.find_first_of(delimiters,1);
	if (string::npos!=pos)
	{
		str=str.substr(0, pos);
	}
}
*/



#endif  // 










