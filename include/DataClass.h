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

#ifndef FQ_DataClass_H_
#define FQ_DataClass_H_

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
//
typedef  long long LLongA;

using namespace std;

bool NArry[256]={false};
bool LowArry[256]={true};
int  n_thread=4;
int  VECMAX =1024*4;

struct OUTIOGZ
{
	ogzstream OUTAPE;
	ogzstream OUTASE;
	ogzstream OUTBPE;
	ogzstream OUTBSE;
	int PE;
	int SEAA;
	int SEBB;
};

struct OUTIO
{
	ofstream OUTAPE;
	ofstream OUTASE;
	ofstream OUTBPE;
	ofstream OUTBSE;
	int PE;
	int SEAA;
	int SEBB;
};

//
class Para_A24 {
	public:
		string InFq1;
		string InFq2;
		string OutFq1;
		string InSeFq ;

//		char N_seq;
		int LowQint ;
		unsigned long ReadNumber;
		int AverQ;
		int MinBaseQ;
		int Kmer;
		int Windows;
		int MinCount;
		int ReadLength;
		int HalfReadLength;
//		int EndReadFlag;
//		int MaxReadPCR ;
//		int EndLastWindows;
		bool OUTGZ ;
		bool PCRA;
		bool FILTER;
		bool OutFa;
		int N_Number;
		double  N_Ration;
		int MinReadKmerCount;
		bool IDAdd;
		bool KmerStatOut;
//		bool PESE;
//		bool allRead;
		Para_A24()
		{
			InFq1="";
			InFq2="";
			OutFq1="";
//			N_seq='N';
			InSeFq="" ;
			LowQint=64;
			MinBaseQ=0;
			AverQ=20;
			Kmer=31;
//			allRead=false;
			Windows=5;
			MinCount=3;
			N_Number=2;
			N_Ration=0.1;
			ReadNumber=1000000;
			ReadLength=0;
			HalfReadLength=0;
//			EndReadFlag=0;
//			MaxReadPCR=4000000;
			//EndLastWindows=60;
			OUTGZ=false;
			PCRA=false;
			OutFa=true;
			MinReadKmerCount=5;
			IDAdd=false;
			KmerStatOut=false;
			FILTER=true;
//			hVal=0;
//			PESE=false;
		}
};

#endif  // 










