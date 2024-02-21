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
#include "gzstream.c"
#include "kseq.h"
//
typedef  long long LLongA;

using namespace std;

bool NArry[256]={false};
bool LowArry[256]={true};
int n_thread=1;
int VECMAX =1024*1024;
int BinWind = VECMAX;
int BATCH_SIZE = BinWind;

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
		int LowQint ;
		int AverQ;
		int MinBaseQ;
		int Kmer;
		int Windows;
		int MinCount;
		int ReadLength;
		int HalfReadLength;
		int PCRA;
		bool OUTGZ ;
		bool FILTER;
		bool OutFa;
		int N_Number;
		int MinReadKmerCount;
		bool KmerStatOut;
		unsigned long ReadNumber;
		double  N_Ration;

		string InFq1;
		string InFq2;
		string OutFq1;
		string InSeFq ;


		Para_A24()
		{
			LowQint=64;
			MinBaseQ=0;
			FILTER=true;
			AverQ=20;
			Kmer=31;
			Windows=5;
			MinCount=3;
			N_Number=2;
			N_Ration=0.1;
			ReadNumber=1000000;
			ReadLength=0;
			HalfReadLength=0;
			PCRA=0;
			OUTGZ=false;
			OutFa=true;
			MinReadKmerCount=5;
			KmerStatOut=false;
			InFq1="";
			InFq2="";
			OutFq1="";
			InSeFq="" ;
		}
};

#endif  // 










