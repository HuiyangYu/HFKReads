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

#ifndef FQ_KmerSM1_H_
#define FQ_KmerSM1_H_

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
#include "./comm.hpp"
//

using namespace std;
//KSEQ_INIT(gzFile, gzread)

int M1_RunFQFilterPE (Para_A24 * P2In)
{
	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  AAASSS ;
	vector <string>  AAAQQQ ;
	vector <string>  AAAIII ;

	vector <string>  BBBSSS ;
	vector <string>  BBBQQQ ;
	vector <string>  BBBIII ;


	int BinWind=VECMAX; 
	int BATCH_SIZE;
	BATCH_SIZE=BinWind*n_thread;

	if (BATCH_SIZE > (P2In->ReadNumber)) 
	{ 
		BinWind =(P2In->ReadNumber)/n_thread; 
		if (BinWind<2) {BinWind=2;} ;
		BATCH_SIZE=BinWind*n_thread;
	}

	AAASSS.resize(BATCH_SIZE+2);
	AAAQQQ.resize(BATCH_SIZE+2);
	AAAIII.resize(BATCH_SIZE+2);

	BBBSSS.resize(BATCH_SIZE+2);
	BBBQQQ.resize(BATCH_SIZE+2);
	BBBIII.resize(BATCH_SIZE+2);

	int A=0;

	std::vector<std::thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];

	bool  *PASS1 =new bool [BATCH_SIZE];
	bool  *PASS2 =new bool [BATCH_SIZE];

	igzstream INA ((P2In->InFq1).c_str(),ifstream::in);
	igzstream INB ((P2In->InFq2).c_str(),ifstream::in);
	INA.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
	INB.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
	string ID_1 ,seq_1,temp_1,Quly_1 ;
	string ID_2 ,seq_2,temp_2,Quly_2 ;

	string OUT=(P2In->OutFq1);

	if (P2In->OUTGZ)
	{
		string outputFileAPE=OUT+"_pe_1.fq.gz";
		string outputFileBPE=OUT+"_pe_2.fq.gz";
		string outputFileASE=OUT+"_se_1.fq.gz";
		string outputFileBSE=OUT+"_se_2.fq.gz";

		if (P2In->OutFa)
		{	
			outputFileAPE=OUT+"_pe_1.fa.gz";
			outputFileBPE=OUT+"_pe_2.fa.gz";
			outputFileASE=OUT+"_se_1.fa.gz";
			outputFileBSE=OUT+"_se_2.fa.gz";
		}

		OUTIOGZ OUTHanDle ;
		OUTHanDle.PE=0;OUTHanDle.SEAA=0;OUTHanDle.PE=0;OUTHanDle.SEBB=0;

		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		OUTHanDle.OUTASE.open(outputFileASE.c_str());
		OUTHanDle.OUTBPE.open(outputFileBPE.c_str());
		OUTHanDle.OUTBSE.open(outputFileBSE.c_str());
		OUTHanDle.OUTAPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTASE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBSE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);

		int CountFQ=0;

		if (P2In->OutFa)
		{
			for (A=0 ; A<(P2In->ReadNumber) && (!INA.eof()) ; A++)
			{
				getline(INA,ID_1);
				getline(INA,seq_1);
				getline(INA,temp_1);
				getline(INA,Quly_1);

				getline(INB,ID_2);
				getline(INB,seq_2);
				getline(INB,temp_2);
				getline(INB,Quly_2);

				if (ID_1.empty())   {continue ; }
				AAASSS[CountFQ]=seq_1;
				AAAQQQ[CountFQ]=Quly_1;
				AAAIII[CountFQ]=ID_1;

				BBBSSS[CountFQ]=seq_2;
				BBBQQQ[CountFQ]=Quly_2;
				BBBIII[CountFQ]=ID_2;

				CountFQ++;
				if (CountFQ==BATCH_SIZE)
				{
					for (int i = 0; i < n_thread; i++)
					{
						Start[i]=i*BinWind;
						End[i]=Start[i]+BinWind;
						threads.push_back(std::thread(M1_FilterFQPE,P2In,PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ),std::ref(BBBSSS),std::ref(BBBQQQ)));
					}

					for (auto& thread : threads)
					{
						thread.join();
					}
					threads.clear();
					RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

					for (int j = 0; j < CountFQ; j++)
					{
						if (PASS1[j]  & PASS2[j] )
						{	
							AAAIII[j][0]='>';
							BBBIII[j][0]='>';
							OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
							OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
							OUTHanDle.PE++;
						}
						else if (PASS1[j])
						{
							AAAIII[j][0]='>';
							OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
							OUTHanDle.SEAA++;
						}
						else if (PASS2[j])
						{
							BBBIII[j][0]='>';
							OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
							OUTHanDle.SEBB++;
						}
					}
					CountFQ=0;
				}
			}

			if (CountFQ!=0)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>CountFQ)  {  End[i]=CountFQ;  }

					threads.push_back(std::thread(M1_FilterFQPE,P2In,PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ),std::ref(BBBSSS),std::ref(BBBQQQ)));
				}
				for (auto& thread : threads)
				{
					thread.join();
				}
				threads.clear();
				RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j]  & PASS2[j] )
					{	
						AAAIII[j][0]='>';
						BBBIII[j][0]='>';
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
						OUTHanDle.PE++;
					}
					else if (PASS1[j])
					{
						AAAIII[j][0]='>';
						OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.SEAA++;
					}
					else if (PASS2[j])
					{
						BBBIII[j][0]='>';
						OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
						OUTHanDle.SEBB++;
					}
				}
				CountFQ=0;
			}

		}
		else
		{

			for (A=0 ; A<(P2In->ReadNumber) && (!INA.eof()) ; A++)
			{
				getline(INA,ID_1);
				getline(INA,seq_1);
				getline(INA,temp_1);
				getline(INA,Quly_1);

				getline(INB,ID_2);
				getline(INB,seq_2);
				getline(INB,temp_2);
				getline(INB,Quly_2);

				if (ID_1.empty())   {continue ; }
				AAASSS[CountFQ]=seq_1;
				AAAQQQ[CountFQ]=Quly_1;
				AAAIII[CountFQ]=ID_1;

				BBBSSS[CountFQ]=seq_2;
				BBBQQQ[CountFQ]=Quly_2;
				BBBIII[CountFQ]=ID_2;

				CountFQ++;
				if (CountFQ==BATCH_SIZE)
				{
					for (int i = 0; i < n_thread; i++)
					{
						Start[i]=i*BinWind;
						End[i]=Start[i]+BinWind;
						threads.push_back(std::thread(M1_FilterFQPE,P2In,PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ),std::ref(BBBSSS),std::ref(BBBQQQ)));
					}

					for (auto& thread : threads)
					{
						thread.join();
					}

					threads.clear();
					RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);


					for (int j = 0; j < CountFQ; j++)
					{
						if (PASS1[j]  & PASS2[j] )
						{						
							OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
							OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
							OUTHanDle.PE++;
						}
						else if (PASS1[j])
						{
							OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
							OUTHanDle.SEAA++;
						}
						else if (PASS2[j])
						{
							OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
							OUTHanDle.SEBB++;
						}
					}
					CountFQ=0;
				}
			}

			if (CountFQ!=0)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>CountFQ)  {  End[i]=CountFQ;  }

					threads.push_back(std::thread(M1_FilterFQPE,P2In,PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ),std::ref(BBBSSS),std::ref(BBBQQQ)));
				}
				for (auto& thread : threads)
				{
					thread.join();
				}
				threads.clear();
				RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j]  & PASS2[j] )
					{						
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
						OUTHanDle.PE++;
					}
					else if (PASS1[j])
					{
						OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OUTHanDle.SEAA++;
					}
					else if (PASS2[j])
					{
						OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
						OUTHanDle.SEBB++;
					}
				}
				CountFQ=0;
			}

		}

		OUTHanDle.OUTAPE.close();
		OUTHanDle.OUTASE.close();
		OUTHanDle.OUTBPE.close();
		OUTHanDle.OUTBSE.close();

		if (!INA.eof()) {  getline(INA,ID_1); if (INA.eof()) { A++;} }
		if (INA.eof()) {A--;	cout<<"INFO: ALL reads "<<A*2<<" are read done"<<endl;	}

		cout<<"INFO: output PE1 read number is "<<OUTHanDle.PE+OUTHanDle.SEAA<<"\n";
		cout<<"INFO: output PE2 read number is "<<OUTHanDle.PE+OUTHanDle.SEBB<<"\n";
		cout<<"INFO: paired PE1 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: paired PE2 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: un-paired PE1 read number is "<<OUTHanDle.SEAA<<"\n";
		cout<<"INFO: un-paired PE2 read number is "<<OUTHanDle.SEBB<<"\n";

	}
	else
	{
		string outputFileAPE=OUT+"_pe_1.fq";
		string outputFileBPE=OUT+"_pe_2.fq";
		string outputFileASE=OUT+"_se_1.fq";
		string outputFileBSE=OUT+"_se_2.fq";

		if (P2In->OutFa)
		{	
			outputFileAPE=OUT+"_pe_1.fa";
			outputFileBPE=OUT+"_pe_2.fa";
			outputFileASE=OUT+"_se_1.fa";
			outputFileBSE=OUT+"_se_2.fa";
		}

		OUTIO OUTHanDle ;
		OUTHanDle.PE=0;OUTHanDle.SEAA=0;OUTHanDle.PE=0;OUTHanDle.SEBB=0;

		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		OUTHanDle.OUTASE.open(outputFileASE.c_str());
		OUTHanDle.OUTBPE.open(outputFileBPE.c_str());
		OUTHanDle.OUTBSE.open(outputFileBSE.c_str());
		OUTHanDle.OUTAPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTASE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBSE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);

		int CountFQ=0;

		if (P2In->OutFa)
		{

			for (A=0 ; A<(P2In->ReadNumber) && (!INA.eof()) ; A++)
			{
				getline(INA,ID_1);
				getline(INA,seq_1);
				getline(INA,temp_1);
				getline(INA,Quly_1);

				getline(INB,ID_2);
				getline(INB,seq_2);
				getline(INB,temp_2);
				getline(INB,Quly_2);

				if (ID_1.empty())   {continue ; }
				AAASSS[CountFQ]=seq_1;
				AAAQQQ[CountFQ]=Quly_1;
				AAAIII[CountFQ]=ID_1;

				BBBSSS[CountFQ]=seq_2;
				BBBQQQ[CountFQ]=Quly_2;
				BBBIII[CountFQ]=ID_2;

				CountFQ++;
				if (CountFQ==BATCH_SIZE)
				{
					for (int i = 0; i < n_thread; i++)
					{
						Start[i]=i*BinWind;
						End[i]=Start[i]+BinWind;
						threads.push_back(std::thread(M1_FilterFQPE,P2In,PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ),std::ref(BBBSSS),std::ref(BBBQQQ)));
					}

					for (auto& thread : threads)
					{
						thread.join();
					}

					threads.clear();
					RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);


					for (int j = 0; j < CountFQ; j++)
					{
						if (PASS1[j]  & PASS2[j] )
						{	
							AAAIII[j][0]='>';
							BBBIII[j][0]='>';
							OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
							OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
							OUTHanDle.PE++;
						}
						else if (PASS1[j])
						{
							AAAIII[j][0]='>';
							OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
							OUTHanDle.SEAA++;
						}
						else if (PASS2[j])
						{
							BBBIII[j][0]='>';
							OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
							OUTHanDle.SEBB++;
						}
					}
					CountFQ=0;
				}
			}

			if (CountFQ!=0)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>CountFQ)  {  End[i]=CountFQ;  }

					threads.push_back(std::thread(M1_FilterFQPE,P2In,PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ),std::ref(BBBSSS),std::ref(BBBQQQ)));
				}
				for (auto& thread : threads)
				{
					thread.join();
				}
				threads.clear();
				RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

				for (int j = 0; j < CountFQ; j++)
				{

					if (PASS1[j]  & PASS2[j] )
					{	
						AAAIII[j][0]='>';
						BBBIII[j][0]='>';
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
						OUTHanDle.PE++;
					}
					else if (PASS1[j])
					{
						AAAIII[j][0]='>';
						OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.SEAA++;
					}
					else if (PASS2[j])
					{
						BBBIII[j][0]='>';
						OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
						OUTHanDle.SEBB++;
					}
				}
				CountFQ=0;
			}

		}
		else
		{

			for (A=0 ; A<(P2In->ReadNumber) && (!INA.eof()) ; A++)
			{
				getline(INA,ID_1);
				getline(INA,seq_1);
				getline(INA,temp_1);
				getline(INA,Quly_1);

				getline(INB,ID_2);
				getline(INB,seq_2);
				getline(INB,temp_2);
				getline(INB,Quly_2);

				if (ID_1.empty())   {continue ; }
				AAASSS[CountFQ]=seq_1;
				AAAQQQ[CountFQ]=Quly_1;
				AAAIII[CountFQ]=ID_1;

				BBBSSS[CountFQ]=seq_2;
				BBBQQQ[CountFQ]=Quly_2;
				BBBIII[CountFQ]=ID_2;

				CountFQ++;
				if (CountFQ==BATCH_SIZE)
				{
					for (int i = 0; i < n_thread; i++)
					{
						Start[i]=i*BinWind;
						End[i]=Start[i]+BinWind;
						threads.push_back(std::thread(M1_FilterFQPE,P2In,PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ),std::ref(BBBSSS),std::ref(BBBQQQ)));
					}

					for (auto& thread : threads)
					{
						thread.join();
					}

					threads.clear();
					RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

					for (int j = 0; j < CountFQ; j++)
					{
						if (PASS1[j]  & PASS2[j] )
						{						
							OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
							OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
							OUTHanDle.PE++;
						}
						else if (PASS1[j])
						{
							OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
							OUTHanDle.SEAA++;
						}
						else if (PASS2[j])
						{
							OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
							OUTHanDle.SEBB++;
						}
					}
					CountFQ=0;
				}
			}

			if (CountFQ!=0)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>CountFQ)  {  End[i]=CountFQ;  }

					threads.push_back(std::thread(M1_FilterFQPE,P2In,PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ),std::ref(BBBSSS),std::ref(BBBQQQ)));
				}
				for (auto& thread : threads)
				{
					thread.join();
				}
				threads.clear();
				RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j]  & PASS2[j] )
					{						
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
						OUTHanDle.PE++;
					}
					else if (PASS1[j])
					{
						OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OUTHanDle.SEAA++;
					}
					else if (PASS2[j])
					{
						OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
						OUTHanDle.SEBB++;
					}
				}
				CountFQ=0;
			}

		}

		OUTHanDle.OUTAPE.close();
		OUTHanDle.OUTASE.close();
		OUTHanDle.OUTBPE.close();
		OUTHanDle.OUTBSE.close();

		if (!INA.eof()) {  getline(INA,ID_1); if (INA.eof()) { A++;} }
		if (INA.eof()) {A--;	cout<<"INFO: ALL reads "<<A*2<<" are read done"<<endl;	}

		cout<<"INFO: output PE1 read number is "<<OUTHanDle.PE+OUTHanDle.SEAA<<"\n";
		cout<<"INFO: output PE2 read number is "<<OUTHanDle.PE+OUTHanDle.SEBB<<"\n";
		cout<<"INFO: paired PE1 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: paired PE2 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: un-paired PE1 read number is "<<OUTHanDle.SEAA<<"\n";
		cout<<"INFO: un-paired PE2 read number is "<<OUTHanDle.SEBB<<"\n";

	}

	INA.close();
	INB.close();
	delete [] Start ;
	delete [] End;
	delete [] PASS1;
	delete [] PASS2;

	return 0;
}

int M1_RunFQFilterSE (Para_A24 * P2In)
{

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  AAASSS ;
	vector <string>  AAAQQQ ;
	vector <string>  AAAIII ;

	int BinWind=VECMAX; 
	int BATCH_SIZE;
	BATCH_SIZE=BinWind*n_thread;

	if (BATCH_SIZE > (P2In->ReadNumber)) 
	{ 
		BinWind =(P2In->ReadNumber)/n_thread; 
		if (BinWind<2) {BinWind=2;} ;
		BATCH_SIZE=BinWind*n_thread;
	}

	AAASSS.resize(BATCH_SIZE+2);
	AAAQQQ.resize(BATCH_SIZE+2);
	AAAIII.resize(BATCH_SIZE+2);

	int A=0;
	int OutPutReadCount=0;

	std::vector<std::thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];

	bool  *PASS =new bool [BATCH_SIZE];

	igzstream INA ((P2In->InFq1).c_str(),ifstream::in);
	INA.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
	string ID_1 ,seq_1,temp_1,Quly_1 ;

	string OUT=(P2In->OutFq1);

	string outputFileAPE=OUT+".fa";

	if (!P2In->OutFa)
	{
		outputFileAPE=OUT+".fq";
	}

	if (P2In->OUTGZ)
	{
		outputFileAPE=outputFileAPE+".gz";			
		OUTIOGZ OUTHanDle ;
		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		int CountFQ=0;

		if (P2In->OutFa)
		{
			for (A=0 ; A<(P2In->ReadNumber) && (!INA.eof()) ; A++)
			{
				getline(INA,ID_1);
				getline(INA,seq_1);
				getline(INA,temp_1);
				getline(INA,Quly_1);

				if (ID_1.empty())   {continue ; }
				AAASSS[CountFQ]=seq_1;
				AAAQQQ[CountFQ]=Quly_1;
				AAAIII[CountFQ]=ID_1;

				CountFQ++;
				if (CountFQ==BATCH_SIZE)
				{
					for (int i = 0; i < n_thread; i++)
					{
						Start[i]=i*BinWind;
						End[i]=Start[i]+BinWind;
						threads.push_back(std::thread(FilterFQSE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ)));
					}

					for (auto& thread : threads)
					{
						thread.join();
					}

					threads.clear();
					RmPCRSE(P2In,PASS,CountFQ,AAASSS);

					for (int j = 0; j < CountFQ; j++)
					{
						if (PASS[j])
						{						
							AAAIII[j][0]='>';
							OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
							OutPutReadCount++;
						}
					}
					CountFQ=0;
				}
			}

			if (CountFQ!=0)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
					threads.push_back(std::thread(FilterFQSE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ)));
				}
				for (auto& thread : threads)
				{
					thread.join();
				}
				threads.clear();
				RmPCRSE(P2In,PASS,CountFQ,AAASSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS[j])
					{
						AAAIII[j][0]='>';
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OutPutReadCount++;
					}
				}
				CountFQ=0;
			}

		}
		else
		{
			for (A=0 ; A<(P2In->ReadNumber) && (!INA.eof()) ; A++)
			{
				getline(INA,ID_1);
				getline(INA,seq_1);
				getline(INA,temp_1);
				getline(INA,Quly_1);

				if (ID_1.empty())   {continue ; }
				AAASSS[CountFQ]=seq_1;
				AAAQQQ[CountFQ]=Quly_1;
				AAAIII[CountFQ]=ID_1;

				CountFQ++;
				if (CountFQ==BATCH_SIZE)
				{
					for (int i = 0; i < n_thread; i++)
					{
						Start[i]=i*BinWind;
						End[i]=Start[i]+BinWind;
						threads.push_back(std::thread(FilterFQSE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ)));
					}

					for (auto& thread : threads)
					{
						thread.join();
					}

					threads.clear();
					RmPCRSE(P2In,PASS,CountFQ,AAASSS);

					for (int j = 0; j < CountFQ; j++)
					{
						if (PASS[j])
						{						
							OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
							OutPutReadCount++;
						}
					}
					CountFQ=0;
				}
			}

			if (CountFQ!=0)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
					threads.push_back(std::thread(FilterFQSE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ)));
				}
				for (auto& thread : threads)
				{
					thread.join();
				}
				threads.clear();
				RmPCRSE(P2In,PASS,CountFQ,AAASSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS[j])
					{
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OutPutReadCount++;
					}
				}
				CountFQ=0;
			}
		}
		OUTHanDle.OUTAPE.close();
	}
	else
	{
		OUTIO OUTHanDle ;
		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());

		int CountFQ=0;

		if (P2In->OutFa)
		{
			for (A=0 ; A<(P2In->ReadNumber) && (!INA.eof()) ; A++)
			{
				getline(INA,ID_1);
				getline(INA,seq_1);
				getline(INA,temp_1);
				getline(INA,Quly_1);

				if (ID_1.empty())   {continue ; }
				AAASSS[CountFQ]=seq_1;
				AAAQQQ[CountFQ]=Quly_1;
				AAAIII[CountFQ]=ID_1;

				CountFQ++;
				if (CountFQ==BATCH_SIZE)
				{
					for (int i = 0; i < n_thread; i++)
					{
						Start[i]=i*BinWind;
						End[i]=Start[i]+BinWind;
						threads.push_back(std::thread(FilterFQSE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ)));
					}

					for (auto& thread : threads)
					{
						thread.join();
					}

					threads.clear();
					RmPCRSE(P2In,PASS,CountFQ,AAASSS);

					for (int j = 0; j < CountFQ; j++)
					{
						if (PASS[j])
						{						
							AAAIII[j][0]='>';
							OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
							OutPutReadCount++;
						}
					}
					CountFQ=0;
				}
			}

			if (CountFQ!=0)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
					threads.push_back(std::thread(FilterFQSE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ)));
				}
				for (auto& thread : threads)
				{
					thread.join();
				}
				threads.clear();
				RmPCRSE(P2In,PASS,CountFQ,AAASSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS[j])
					{
						AAAIII[j][0]='>';
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OutPutReadCount++;
					}
				}
				CountFQ=0;
			}

		}
		else
		{
			for (A=0 ; A<(P2In->ReadNumber) && (!INA.eof()) ; A++)
			{
				getline(INA,ID_1);
				getline(INA,seq_1);
				getline(INA,temp_1);
				getline(INA,Quly_1);

				if (ID_1.empty())   {continue ; }
				AAASSS[CountFQ]=seq_1;
				AAAQQQ[CountFQ]=Quly_1;
				AAAIII[CountFQ]=ID_1;

				CountFQ++;
				if (CountFQ==BATCH_SIZE)
				{
					for (int i = 0; i < n_thread; i++)
					{
						Start[i]=i*BinWind;
						End[i]=Start[i]+BinWind;
						threads.push_back(std::thread(FilterFQSE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ)));
					}

					for (auto& thread : threads)
					{
						thread.join();
					}

					threads.clear();
					RmPCRSE(P2In,PASS,CountFQ,AAASSS);

					for (int j = 0; j < CountFQ; j++)
					{
						if (PASS[j])
						{						
							OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
							OutPutReadCount++;
						}
					}
					CountFQ=0;
				}
			}

			if (CountFQ!=0)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
					threads.push_back(std::thread(FilterFQSE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ)));
				}
				for (auto& thread : threads)
				{
					thread.join();
				}
				threads.clear();
				RmPCRSE(P2In,PASS,CountFQ,AAASSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS[j])
					{
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OutPutReadCount++;
					}
				}
				CountFQ=0;
			}
		}
		OUTHanDle.OUTAPE.close();
	}

	if (!INA.eof()) {  getline(INA,ID_1); if (INA.eof()) { A++;} }
	if (INA.eof()) {A--;	cout<<"INFO: ALL reads "<<A<<" are read done"<<endl;	}

	INA.close();

	delete [] Start ;
	delete [] End;
	delete [] PASS;
	cout<<"INFO: output SE read number is "<<OutPutReadCount<<endl;
	return 0;
}

int M1_RunFAFilterPE (Para_A24 * P2In)
{
	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  AAASSS ;
	vector <string>  AAAIII ;

	vector <string>  BBBSSS ;
	vector <string>  BBBIII ;

	int BinWind=VECMAX; 
	int BATCH_SIZE;
	BATCH_SIZE=BinWind*n_thread;

	if (BATCH_SIZE > (P2In->ReadNumber)) 
	{ 
		BinWind =(P2In->ReadNumber)/n_thread; 
		if (BinWind<2) {BinWind=2;} ;
		BATCH_SIZE=BinWind*n_thread;
	}

	AAASSS.resize(BATCH_SIZE+2);
	AAAIII.resize(BATCH_SIZE+2);

	BBBSSS.resize(BATCH_SIZE+2);
	BBBIII.resize(BATCH_SIZE+2);

	std::vector<std::thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];

	bool  *PASS1 =new bool [BATCH_SIZE];
	bool  *PASS2 =new bool [BATCH_SIZE];

	gzFile fpAA;
	kseq_t *seqAA;

	gzFile fpBB;
	kseq_t *seqBB;

 	fpAA = gzopen((P2In->InFq1).c_str(), "r");
	seqAA = kseq_init(fpAA);

	fpBB = gzopen((P2In->InFq2).c_str(), "r");
	seqBB = kseq_init(fpBB);

	int AA ; int A=0;
	string ID_1,ID_2;
	string seqStrAA,seqStrBB;

	string OUT=(P2In->OutFq1);

	string outputFileAPE=OUT+"_pe_1.fa";
	string outputFileBPE=OUT+"_pe_2.fa";
	string outputFileASE=OUT+"_se_1.fa";
	string outputFileBSE=OUT+"_se_2.fa";

	if (P2In->OUTGZ)
	{
		outputFileAPE=OUT+"_pe_1.fa.gz";
		outputFileBPE=OUT+"_pe_2.fa.gz";
		outputFileASE=OUT+"_se_1.fa.gz";
		outputFileBSE=OUT+"_se_2.fa.gz";
		
		OUTIOGZ OUTHanDle ;
		OUTHanDle.PE=0;OUTHanDle.SEAA=0;OUTHanDle.PE=0;OUTHanDle.SEBB=0;
		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		OUTHanDle.OUTASE.open(outputFileASE.c_str());
		OUTHanDle.OUTBPE.open(outputFileBPE.c_str());
		OUTHanDle.OUTBSE.open(outputFileBSE.c_str());
		OUTHanDle.OUTAPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTASE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBSE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);

		int CountFQ=0;

		for (A=0 ; A<(P2In->ReadNumber) && ( (AA = kseq_read(seqAA)) >= 0)  ; A++)
		{

			ID_1=(seqAA->name.s);
			seqStrAA=(seqAA->seq.s);

			AAASSS[CountFQ]=seqStrAA;
			AAAIII[CountFQ]=ID_1;

			AA = kseq_read(seqBB);

			ID_2=(seqBB->name.s);
			seqStrBB=(seqBB->seq.s);
			BBBSSS[CountFQ]=seqStrBB;
			BBBIII[CountFQ]=ID_2;
			CountFQ++;

			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(M1_FilterFAPE,P2In,PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j]  & PASS2[j] )
					{	
						AAAIII[j][0]='>';
						BBBIII[j][0]='>';
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
						OUTHanDle.PE++;
					}
					else if (PASS1[j])
					{
						AAAIII[j][0]='>';
						OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.SEAA++;
					}
					else if (PASS2[j])
					{
						BBBIII[j][0]='>';
						OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
						OUTHanDle.SEBB++;
					}

				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(M1_FilterFAPE,P2In,PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS1[j]  & PASS2[j] )
				{	
					AAAIII[j][0]='>';
					BBBIII[j][0]='>';
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
					OUTHanDle.PE++;
				}
				else if (PASS1[j])
				{
					AAAIII[j][0]='>';
					OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OUTHanDle.SEAA++;
				}
				else if (PASS2[j])
				{
					BBBIII[j][0]='>';
					OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
					OUTHanDle.SEBB++;
				}

			}
			CountFQ=0;
		}

		OUTHanDle.OUTAPE.close();
		OUTHanDle.OUTASE.close();
		OUTHanDle.OUTBPE.close();
		OUTHanDle.OUTBSE.close();

		if  (A==(P2In->ReadNumber) )
		{
			if  ((AA = kseq_read(seqAA)) >= 0) {  }
			else
			{
				cout<<"INFO: ALL reads "<<A*2<<" are read done"<<endl;
			}
		}
		else
		{
			cout<<"INFO: ALL reads "<<A*2<<" are read done"<<endl;
		}

		cout<<"INFO: output PE1 read number is "<<OUTHanDle.PE+OUTHanDle.SEAA<<"\n";
		cout<<"INFO: output PE2 read number is "<<OUTHanDle.PE+OUTHanDle.SEBB<<"\n";
		cout<<"INFO: paired PE1 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: paired PE2 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: un-paired PE1 read number is "<<OUTHanDle.SEAA<<"\n";
		cout<<"INFO: un-paired PE2 read number is "<<OUTHanDle.SEBB<<"\n";

	}
	else
	{
		OUTIO OUTHanDle ;
		OUTHanDle.PE=0;OUTHanDle.SEAA=0;OUTHanDle.PE=0;OUTHanDle.SEBB=0;
		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		OUTHanDle.OUTASE.open(outputFileASE.c_str());
		OUTHanDle.OUTBPE.open(outputFileBPE.c_str());
		OUTHanDle.OUTBSE.open(outputFileBSE.c_str());
		OUTHanDle.OUTAPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTASE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBSE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);

		int CountFQ=0;

		for (A=0 ; A<(P2In->ReadNumber) && ( (AA = kseq_read(seqAA)) >= 0)  ; A++)
		{

			ID_1=(seqAA->name.s);
			seqStrAA=(seqAA->seq.s);

			AAASSS[CountFQ]=seqStrAA;
			AAAIII[CountFQ]=ID_1;

			AA = kseq_read(seqBB);

			ID_2=(seqBB->name.s);
			seqStrBB=(seqBB->seq.s);
			BBBSSS[CountFQ]=seqStrBB;
			BBBIII[CountFQ]=ID_2;
			CountFQ++;

			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(M1_FilterFAPE,P2In,PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j]  & PASS2[j] )
					{	
						AAAIII[j][0]='>';
						BBBIII[j][0]='>';
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
						OUTHanDle.PE++;
					}
					else if (PASS1[j])
					{
						AAAIII[j][0]='>';
						OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.SEAA++;
					}
					else if (PASS2[j])
					{
						BBBIII[j][0]='>';
						OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
						OUTHanDle.SEBB++;
					}

				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(M1_FilterFAPE,P2In,PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS1[j]  & PASS2[j] )
				{	
					AAAIII[j][0]='>';
					BBBIII[j][0]='>';
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
					OUTHanDle.PE++;
				}
				else if (PASS1[j])
				{
					AAAIII[j][0]='>';
					OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OUTHanDle.SEAA++;
				}
				else if (PASS2[j])
				{
					BBBIII[j][0]='>';
					OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
					OUTHanDle.SEBB++;
				}

			}
			CountFQ=0;
		}

		OUTHanDle.OUTAPE.close();
		OUTHanDle.OUTASE.close();
		OUTHanDle.OUTBPE.close();
		OUTHanDle.OUTBSE.close();

		if  (A==(P2In->ReadNumber) )
		{
			if  ((AA = kseq_read(seqAA)) >= 0) {  }
			else
			{
				cout<<"INFO: ALL reads "<<A*2<<" are read done"<<endl;
			}
		}
		else
		{
			cout<<"INFO: ALL reads "<<A*2<<" are read done"<<endl;
		}

		cout<<"INFO: output PE1 read number is "<<OUTHanDle.PE+OUTHanDle.SEAA<<"\n";
		cout<<"INFO: output PE2 read number is "<<OUTHanDle.PE+OUTHanDle.SEBB<<"\n";
		cout<<"INFO: paired PE1 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: paired PE2 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: un-paired PE1 read number is "<<OUTHanDle.SEAA<<"\n";
		cout<<"INFO: un-paired PE2 read number is "<<OUTHanDle.SEBB<<"\n";

	}

	kseq_destroy(seqAA);
	gzclose(fpAA);

	kseq_destroy(seqBB);
	gzclose(fpBB);

	delete [] Start ;
	delete [] End;
	delete [] PASS1;
	delete [] PASS2;

	return 0;
}

int M1_RunFAFilterSE (Para_A24 * P2In)
{

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  AAASSS ;
	vector <string>  AAAIII ;

	int BinWind=VECMAX; 
	int BATCH_SIZE;
	BATCH_SIZE=BinWind*n_thread;

	if (BATCH_SIZE > (P2In->ReadNumber)) 
	{ 
		BinWind =(P2In->ReadNumber)/n_thread; 
		if (BinWind<2) {BinWind=2;} ;
		BATCH_SIZE=BinWind*n_thread;
	}
	
	AAASSS.resize(BATCH_SIZE+2);
	AAAIII.resize(BATCH_SIZE+2);

	std::vector<std::thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];

	bool  *PASS =new bool [BATCH_SIZE];

	gzFile fpAA;
	kseq_t *seqAA;

 	fpAA = gzopen((P2In->InFq1).c_str(), "r");
	seqAA = kseq_init(fpAA);

	int AA ; int A=0;
	int OutPutReadCount=0;
	string ID_1;
	string seqStrAA;

	string OUT=(P2In->OutFq1);
	string outputFileAPE=OUT+".fa";

	if (P2In->OUTGZ)
	{
		outputFileAPE=outputFileAPE+".gz";
		OUTIOGZ OUTHanDle ;
		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		int CountFQ=0;

		for (A=0 ; A<(P2In->ReadNumber) && ( (AA = kseq_read(seqAA)) >= 0)  ; A++)
		{

			ID_1=(seqAA->name.s);
			seqStrAA=(seqAA->seq.s);

			AAASSS[CountFQ]=seqStrAA;
			AAAIII[CountFQ]=ID_1;

			CountFQ++;

			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterFASE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRSE(P2In,PASS,CountFQ,AAASSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS[j])
					{						
						AAAIII[j][0]='>';
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OutPutReadCount++;
					}
				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(FilterFASE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRSE(P2In,PASS,CountFQ,AAASSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS[j])
				{
					AAAIII[j][0]='>';
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OutPutReadCount++;
				}
			}
			CountFQ=0;
		}

		OUTHanDle.OUTAPE.close();

		if  (A==(P2In->ReadNumber) )
		{
			if  ((AA = kseq_read(seqAA)) >= 0) {  }
			else
			{
				cout<<"INFO: ALL reads "<<A<<" are read done"<<endl;
			}
		}
		else
		{
			cout<<"INFO: ALL reads "<<A<<" are read done"<<endl;
		}
	}
	else
	{
		OUTIO OUTHanDle ;
		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		int CountFQ=0;

		for (A=0 ; A<(P2In->ReadNumber) && ( (AA = kseq_read(seqAA)) >= 0)  ; A++)
		{

			ID_1=(seqAA->name.s);
			seqStrAA=(seqAA->seq.s);

			AAASSS[CountFQ]=seqStrAA;
			AAAIII[CountFQ]=ID_1;

			CountFQ++;

			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterFASE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRSE(P2In,PASS,CountFQ,AAASSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS[j])
					{						
						AAAIII[j][0]='>';
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OutPutReadCount++;
					}
				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(FilterFASE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRSE(P2In,PASS,CountFQ,AAASSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS[j])
				{
					AAAIII[j][0]='>';
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OutPutReadCount++;
				}
			}
			CountFQ=0;
		}

		OUTHanDle.OUTAPE.close();

		if  (A==(P2In->ReadNumber) )
		{
			if  ((AA = kseq_read(seqAA)) >= 0) {  }
			else
			{
				cout<<"INFO: ALL reads "<<A<<" are read done"<<endl;
			}
		}
		else
		{
			cout<<"INFO: ALL reads "<<A<<" are read done"<<endl;
		}

	}

	kseq_destroy(seqAA);
	gzclose(fpAA);

	delete [] Start ;
	delete [] End;
	delete [] PASS;
	cout<<"INFO: output SE read number is "<<OutPutReadCount<<endl;

	return 0;
}

//////////////////RunM1_Work///////////////////
int RunM1_Work (Para_A24 * P2In ,int & InPESE )
{
	if (P2In->KmerStatOut)
	{
		cout<<"Warning: Para -f do not work with the [ -m 1 ]\n";
	}

	if ((InPESE==1) && ((P2In->LowQint)!=0) )
	{
		cout <<"INFO: Run paired-end fastq (-m 1) with Phred is "<<(P2In->LowQint)<<endl;
		P2In->ReadNumber=(P2In->ReadNumber)/2;
		M1_RunFQFilterPE(P2In);
	}
	else if ((InPESE==2) && ((P2In->LowQint)!=0) )
	{
		cout <<"INFO: Run single-end fastq (-m 1) with Phred is "<<(P2In->LowQint)<<endl;
		M1_RunFQFilterSE(P2In);
	}
	else if ((InPESE==2) && ((P2In->LowQint)==0) )
	{
		cout <<"INFO: Run single-end fasta (-m 1 )"<<endl;
		M1_RunFAFilterSE (P2In);
	}
	else if ((InPESE==1) && ((P2In->LowQint)==0) )
	{
		cout <<"INFO: Run paired-end fasta (-m 1)"<<endl;
		P2In->ReadNumber=(P2In->ReadNumber)/2;
		M1_RunFAFilterPE(P2In);
	}

	return  0;
}
#endif  //
