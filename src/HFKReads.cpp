/* The GPL-3.0 License

   Copyright (c) 2023- by Huiyang Yu, Weiming He, Chunmei Shi.

*/

#ifndef FQ_KmerSplit_H_
#define FQ_KmerSplit_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <thread>
#include <algorithm>
#include <atomic>
#include <unordered_set>
#include <cmath>
#include <cstdio>
#include <vector>
#include <zlib.h>
#include "kseq.h"
#include "kc-c4-window.c"

typedef  long long LLongA;

using namespace std;

//KSEQ_INIT(gzFile, gzread)

int  print_usage_HFKreads() {
	cout <<""
		"Usage: hfkreads -1 PE1.fq.gz -2 PE2.fq.gz -o OutPrefix\n"
		" Input/Output options:\n"
		"   -1	<str>   paired-end fasta/q file1\n"
		"   -2	<str>   paired-end fasta/q file2\n"
		"   -s	<str>   single-end fasta/q file\n"
		"   -o	<str>   prefix of output file\n"
		" Filter options:\n"
		"   -b	<int>   min base quality [0]\n"
		"   -q	<int>   min average base quality [20]\n"
		"   -l	<int>   min length of read [half]\n"
		"   -r	<float> max unknown base (N) ratio [0.1]\n"
		"   -k	<int>   kmer length [31]\n"
		"   -w	<int>   window size [10]\n"
		"   -m	<int>   min count for high frequency kmer (HFK) [3] \n"
		"   -x	<float>	min ratio of HFK in the read [0.9]\n"
		"   -n	<int>   read number to use [1000000]\n"
		"   -a	        use all the read\n"
		" Other options:\n"
		"   -d           drop the duplicated reads/pairs\n"
		"   -f           output the kmer frequency file\n"
		"   -A           keep base quality in output\n"
		//"   -c           compress the output File\n"
		"   -t           number of threads [1]\n"
		"   -h           show help [v2.03]\n"
		"\n";
	return 1;
}
//

bool NArry[256]={false};
bool LowArry[256]={true};
int n_thread=1;
int VECMAX =102400;
int BinWind = VECMAX;
int BATCH_SIZE = BinWind;

class Para_A24 {
	public:
		int MaxQV;
		int AverQ;
		int MinBaseQ;
		int Kmer;
		int Windows;
		int MinCount;
		int ReadLength;
		int MinLen;
		bool Dedup;
		bool OUTGZ ;
		bool FILTER_LQ;
		bool OutFa;
		int N_Number;
		double MinReadKmerRatio;
		bool KmerStatOut;
		unsigned long ReadNumber;
		double  N_Ratio;

		string InFq1;
		string InFq2;
		string InSeFq ;
		string OutPrefix;
		
		Para_A24() {
			MaxQV=64;
			MinBaseQ=0;
			AverQ=20;
			Kmer=31;
			Windows=10;
			MinCount=3;
			N_Number=2;
			N_Ratio=0.1;
			MinReadKmerRatio=0.9;
			ReadNumber=1000000;
			ReadLength=0;
			MinLen=0;
			Dedup=false;
			OUTGZ=false;
			OutFa=true;
			FILTER_LQ=true;
			KmerStatOut=false;
			InFq1="";
			InFq2="";
			InSeFq="" ;
			OutPrefix="";
		}
};

inline void  LogLackArg(string flag) {
	cerr << "Error: Lack Argument for [ -"<<flag<<" ]"<<endl;
}

string & replace_all(string & str, const string & pattern, const string & replacement) {
	while(true) {
		string::size_type  pos(0);
		if((pos=str.find(pattern))!=string::npos){
			str.replace(pos,pattern.length(),replacement);
		} else{
			break;
		}
	}
	return str;
}

int parse_cmd_HFKreads(int argc, char **argv, Para_A24 * P2In) {
	if (argc <=3) {print_usage_HFKreads();return 0;}

	int err_flag = 0;

	for(int i = 1; i < argc || err_flag; i++) {
		if(argv[i][0] != '-') {
			cerr << "Error: command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		//Input/Output options
		if (flag == "1" ) {
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			P2In->InFq1=argv[i];
		}
		else if (flag == "2") {
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			P2In->InFq2=argv[i];
		}
		else if (flag == "s" ) {
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			P2In->InSeFq=argv[i];
		}
		else if (flag == "o" ) {
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			P2In->OutPrefix=argv[i];
		}

		//Filter low quality
		else if (flag == "b") {
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			P2In->MinBaseQ=atoi(argv[i]);
		}
		else if (flag == "q") {
			if(i + 1 == argc) { LogLackArg(flag) ; return 0;}
			i++;
			P2In->AverQ=atoi(argv[i]);
		}
		else if (flag == "l") {
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->MinLen=atoi(argv[i]);
			if ((P2In->MinLen)<11) {
				P2In->MinLen=11;
				cerr<<"Warings: -l should >= 11, we set it to 11\n";
			}
		}
		else if (flag == "r") {
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->N_Ratio=atof(argv[i]);
		}

		//Filter low freq reads
		else if (flag == "k") {
			if(i + 1 == argc) { LogLackArg(flag) ; return 0;}
			i++;
			P2In->Kmer=atoi(argv[i]);
			if ( (P2In->Kmer) >31 || (P2In->Kmer) < 9) {
				cerr<<"Warning -k  [9,31] ,we modify it to be  31\n ";
				P2In->Kmer=31;
			}
		}
		else if (flag == "w") {
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			P2In->Windows=atoi(argv[i]);
			if  ((P2In->Windows)<1){
				P2In->Windows=1;
				cerr<<"Warnings: -w must be >=1, we change it to 1 \n";
			}
		}
		else if (flag == "m"){
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			P2In->MinCount=atoi(argv[i]);
			if  ((P2In->MinCount)<1) {
				cerr<<"Warnings: -m should >= 1, we set it to 1\n";
				(P2In->MinCount)=1;
			}
		}
		else if (flag == "x") {
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			P2In->MinReadKmerRatio=atof(argv[i]);
			if( (P2In->MinReadKmerRatio)>1 ) {P2In->MinReadKmerRatio=1;}
		}
		else if (flag == "n") {
			if(i + 1 == argc) { LogLackArg(flag) ; return 0;}
			i++;
			P2In->ReadNumber=atoi(argv[i]);
		}
		else if (flag == "A") {
			P2In->OutFa=false;
		}

		//other options
		else if (flag  == "d") {
			P2In->Dedup=true;
		}
		else if (flag  == "c") {
			P2In->OUTGZ=true;
		}
		else if (flag  == "f") {
			P2In->KmerStatOut=true;
		}
		else if (flag  == "a") {
			P2In->ReadNumber=LONG_MAX;
			//P2In->allRead=true;
		}
		else if (flag  ==  "t") {
			if(i + 1 == argc) {LogLackArg( flag ) ; return 0;}
			i++;
			n_thread=atoi(argv[i]);
		}
		else if (flag  == "help" || flag  == "h") {
			print_usage_HFKreads();return 0;
		}
		//
		else if (flag  ==  "u") {
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			VECMAX=atoi(argv[i]);
		}
		else {
			cerr << "Error: UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	//check file type; 0 for unknow; 1 for PE; 2 for SE 
	if ((P2In->InSeFq).empty()){
		if ((P2In->InFq1).empty() || (P2In->InFq2).empty()){
			cerr<< "Error: PE reads should input together"<<endl;
			return 0;
		}else{
			if (access((P2In->InFq1).c_str(), 0) != 0){
				cerr<<"Error: can't find this file -1 "<<(P2In->InFq1)<<endl;
				return 0 ;
			}else if ((access((P2In->InFq2).c_str(), 0) != 0)){
				cerr<<"Error: can't find this file -2 "<<(P2In->InFq2)<<endl;
				return 0 ;
			}else{
				return 1;
			}
		}
	}else{
		if (!(P2In->InFq1).empty() || !(P2In->InFq2).empty()){
			cerr<< "Error: PE reads can't input with SE reads"<<endl;
			return 0;
		}else{
			if (access((P2In->InSeFq).c_str(), 0) != 0){
				cerr<<"Error: can't find this file -s "<<(P2In->InSeFq)<<endl;
				return 0 ;
			}else{
				return 2;
			}
		}
	}
}

int Get_qType(string FilePath, Para_A24 * P2In){

	gzFile fp;
  	kseq_t *seq;
  	int len;

	fp = gzopen((FilePath).c_str(), "r");
  	seq = kseq_init(fp);

	int seqNum=0;
	int maxSeq=5000;
	int minQ=50000;
	int maxQ=0;
	int Lengths[maxSeq];
	string qual;
	for (int A=0 ; A<(maxSeq) && ((len = kseq_read(seq)) >= 0); A++){
		Lengths[seqNum]=(seq->seq.l);
		seqNum++;
		qual=seq->qual.s;
		for (int i=0; i<(seq->qual.l); i++) {
			if(minQ>qual[i]) {
				minQ=qual[i];
			}
			if(maxQ<qual[i]) {
				maxQ=qual[i];
			}
		}
	}
	kseq_destroy(seq);
  	gzclose(fp);

	sort(Lengths, Lengths + seqNum);
	int middleIndex = seqNum / 2;
	P2In->ReadLength = Lengths[middleIndex];
	P2In->N_Number=int((P2In->ReadLength)*(P2In->N_Ratio));

	if ((P2In->MinLen) < ((P2In->Kmer)+1)) {
		P2In->MinLen=(P2In->ReadLength)/2;
		if ((P2In->MinLen)<((P2In->Kmer)+1)){
			(P2In->MinLen)=(P2In->Kmer)+1;
		}
	}
	
	int qType=0;
	if (maxQ>0){
		if(minQ >= 33 &&  minQ <= 78  &&  maxQ >= 33 && maxQ <= 78) {
			qType=33;
		}
		else if (minQ >= 64  &&  minQ <= 108  &&  maxQ >= 64 && maxQ <= 108){
			qType=64;
		}
		else if (minQ < 55) {
			qType=33;
		}
		else {
			qType=64;
		}
		P2In->AverQ=(P2In->AverQ)+qType;
		P2In->MinBaseQ=(P2In->MinBaseQ)+qType;
	}

	return maxQ;
}

inline bool Filter_fq_reads(Para_A24 * P2In, string & seq, string & qual){
	int sumQ=0;
	bool Low=true;
	int NN=0;
	int read_length=seq.length();

	if (read_length< (P2In->MinLen)) {
		return false;
	}

	for(int ix=0 ; ix<read_length ; ix++){
		NN+=NArry[seq[ix]] ;	
		Low=Low & LowArry[qual[ix]];
		sumQ+=qual[ix];
	}

	if (NN>(P2In->N_Number)) {
		return false;
	}
	if (!Low) {
		return false;
	}

	if ((sumQ/read_length) < (P2In->AverQ)) {
		return false;
	}

	return true;
}

inline bool Filter_fa_reads(Para_A24 * P2In, string & seq) {
	int NN=0;
	int read_length=seq.length();
	if (read_length< (P2In->MinLen)) {return false;}
	for(int ix=0; ix<read_length ; ix++) {
		NN+=NArry[seq[ix]];	
	}
	if (NN>(P2In->N_Number)) {return false ;}
	return true ;
}

void Filter_SE_low_qual_reads(Para_A24 * P2In, bool * PASS, int & Start, int & End, 
							  vector <string> & SEQ, vector <string> & QUAL){
	if(QUAL[0].length()<=2){
		for (int i=Start; i<End; i++){
			PASS[i]=Filter_fa_reads(P2In, SEQ[i]);
		}
	} else {
		for (int i=Start; i<End; i++){
			PASS[i]=Filter_fq_reads(P2In, SEQ[i], QUAL[i]);
		}
	}

}

void Filter_PE_low_qual_reads(Para_A24 * P2In, bool * PASS, int & Start, int & End, 
							  vector <string> & PE1_SEQ, vector <string> & PE1_QUAL,
							  vector <string> & PE2_SEQ, vector <string> & PE2_QUAL){
	if(PE1_QUAL[0].length()<=2){
		for (int i=Start; i<End; i++){
			if (Filter_fa_reads(P2In, PE1_SEQ[i]) && Filter_fa_reads(P2In, PE2_SEQ[i])){
				PASS[i]=true;
			} else{
				PASS[i]=false;
			}
		}
	} else {
		for (int i=Start; i<End; i++){
			if(Filter_fq_reads(P2In, PE1_SEQ[i], PE1_QUAL[i]) && Filter_fq_reads(P2In, PE2_SEQ[i], PE2_QUAL[i])){
				PASS[i]=true;
			} else{
				PASS[i]=false;
			}
		}
	}
}

void RmPCR_PE(Para_A24 * P2In, bool * PASS, int & End, vector <string> & PE1_SEQ, vector <string> & PE2_SEQ)
{
	if (P2In->Dedup) {
		string Cat;
		unordered_set<string> seen_seqs;
		for (int ii = 0; ii < End; ii++) {
			if (PASS[ii]) {
				Cat = PE1_SEQ[ii] + PE2_SEQ[ii];
				if (seen_seqs.find(Cat) != seen_seqs.end()) {
					PASS[ii] = false;
				} else {
					seen_seqs.insert(Cat);
				}
			}
		}
	}
}

void RmPCR_SE(Para_A24 * P2In, bool * PASS, int & End, vector <string> & SEQ)
{
    if(P2In->Dedup){
		unordered_set<string> seen_seqs;
        for (int ii=0; ii<End; ii++) {
            if (PASS[ii]) {
                if (seen_seqs.find(SEQ[ii]) != seen_seqs.end()) {
                    PASS[ii] = false;
                } else {
                    seen_seqs.insert(SEQ[ii]);
                }
            }
        }
    }
}

int Out_SE_seq(Para_A24 * P2In, int &seq_num, ofstream &OUTH, bool *PASS, vector <string> &ID, vector <string> &SEQ, vector <string> &QUAL){
	int out_number=0;
	if (P2In->OutFa){
		for (int j = 0; j < seq_num; j++) {
			if (PASS[j]) {						
				OUTH << ">"<<ID[j]<<"\n"<<SEQ[j]<<"\n";
				out_number++;
			}
		}
	} else {
		for (int j = 0; j < seq_num; j++) {
			if (PASS[j]) {
				OUTH<< "@"<<ID[j]<<"\n"<<SEQ[j]<<"\n+\n"<<QUAL[j]<<"\n";
				out_number++;
			}
		}
	}
	return out_number;
}

vector<int> Out_PE_seq(Para_A24 * P2In, int &seq_num, ofstream &OUTH_PE1, ofstream &OUTH_PE2,
			  ofstream &OUTH_SE1, ofstream &OUTH_SE2, bool *PASS_PE1, bool *PASS_PE2,
			  vector <string> &PE1_ID, vector <string> &PE1_SEQ, vector <string> &PE1_QUAL,
			  vector <string> &PE2_ID, vector <string> &PE2_SEQ, vector <string> &PE2_QUAL){
	int pe_number=0;
	int se_number=0;
	if (P2In->OutFa){
		for (int j = 0; j < seq_num; j++) {
			if (PASS_PE1[j] && PASS_PE2[j]) {						
				OUTH_PE1 << ">"<< PE1_ID[j]<<"\n"<<PE1_SEQ[j]<<"\n";
				OUTH_PE2 << ">"<< PE2_ID[j]<<"\n"<<PE2_SEQ[j]<<"\n";
				pe_number++;
			} else if(PASS_PE1[j]){
				OUTH_SE1 << ">"<< PE1_ID[j]<<"\n"<<PE1_SEQ[j]<<"\n";
				se_number++;
			} else if(PASS_PE2[j]){
				OUTH_SE2 << ">"<< PE2_ID[j]<<"\n"<<PE2_SEQ[j]<<"\n";
			}
		}
	} else {
		for (int j = 0; j < seq_num; j++) {
			if (PASS_PE1[j] && PASS_PE2[j]) {
				OUTH_PE1<< "@"<<PE1_ID[j]<<"\n"<<PE1_SEQ[j]<<"\n+\n"<<PE1_QUAL[j]<<"\n";
				OUTH_PE2<< "@"<<PE2_ID[j]<<"\n"<<PE2_SEQ[j]<<"\n+\n"<<PE2_QUAL[j]<<"\n";
				pe_number++;
			}else if(PASS_PE1[j]){
				OUTH_SE1<< "@"<<PE1_ID[j]<<"\n"<<PE1_SEQ[j]<<"\n+\n"<<PE1_QUAL[j]<<"\n";
				se_number++;
			}else if(PASS_PE2[j]){
				OUTH_SE2<< "@"<<PE2_ID[j]<<"\n"<<PE2_SEQ[j]<<"\n+\n"<<PE2_QUAL[j]<<"\n";
			}
		}
	}
	vector<int> out_number = {pe_number, se_number};
	return out_number;
}

int Run_PE_low_qual_batch(Para_A24 * P2In, int &seq_num, ofstream &OUTH_PE1, ofstream &OUTH_PE2,
				vector <string> &PE1_ID, vector <string> &PE1_SEQ, vector <string> &PE1_QUAL,
				vector <string> &PE2_ID, vector <string> &PE2_SEQ, vector <string> &PE2_QUAL){
	
	std::vector<std::thread> threads;
	int * Start =new int [n_thread];
	int * End =new int [n_thread];
	bool  *PASS =new bool [BATCH_SIZE];

	for (int i = 0; i < n_thread; i++){
		Start[i]=i*BinWind;
		End[i]=Start[i]+BinWind;
		if (End[i]>seq_num) {End[i]=seq_num;} if (Start[i]>=End[i]) {continue;}
		threads.push_back(thread(Filter_PE_low_qual_reads,P2In,PASS,ref(Start[i]),ref(End[i]),
						         ref(PE1_SEQ),ref(PE1_QUAL),ref(PE2_SEQ),ref(PE2_QUAL)));
	}

	for (auto& thread : threads) {
		thread.join();
	}

	threads.clear();
	RmPCR_PE(P2In,PASS,seq_num,PE1_SEQ,PE2_SEQ);
	int out_number=Out_SE_seq(P2In,seq_num,OUTH_PE1,PASS,PE1_ID,PE1_SEQ,PE1_QUAL);
	Out_SE_seq(P2In,seq_num,OUTH_PE2,PASS,PE2_ID,PE2_SEQ,PE2_QUAL);

	delete [] Start;
	delete [] End;
	delete [] PASS;

	return out_number;
}

int Run_SE_low_qual_batch(Para_A24 * P2In, int &seq_num, ofstream &OUTH_SE, 
				vector <string> &SE_ID, vector <string> &SE_SEQ, vector <string> &SE_QUAL){
	
	std::vector<std::thread> threads;
	int * Start =new int [n_thread];
	int * End =new int [n_thread];
	bool  *PASS =new bool [BATCH_SIZE];

	for (int i = 0; i < n_thread; i++){
		Start[i]=i*BinWind;
		End[i]=Start[i]+BinWind;
		if (End[i]>seq_num) {End[i]=seq_num;} if (Start[i]>=End[i]) {continue;}
		threads.push_back(thread(Filter_SE_low_qual_reads,P2In,PASS,
								 ref(Start[i]),ref(End[i]), ref(SE_SEQ),ref(SE_QUAL)));
	}

	for (auto& thread : threads) {
		thread.join();
	}

	threads.clear();

	RmPCR_SE(P2In, PASS, seq_num, SE_SEQ);
	int out_number=Out_SE_seq(P2In,seq_num,OUTH_SE,PASS,SE_ID,SE_SEQ,SE_QUAL);

	delete [] Start;
	delete [] End;
	delete [] PASS;

	return out_number;
}

int Run_PE_low_qual_filter (Para_A24 * P2In, vector<std::string> & FilePath) {

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  PE1_ID;
	vector <string>  PE1_SEQ;
	vector <string>  PE1_QUAL;
	
	vector <string>  PE2_ID ;
	vector <string>  PE2_SEQ ;
	vector <string>  PE2_QUAL ;

	PE1_ID.resize(BATCH_SIZE+2);
	PE1_SEQ.resize(BATCH_SIZE+2);
	PE1_QUAL.resize(BATCH_SIZE+2);
	
	PE2_ID.resize(BATCH_SIZE+2);
	PE2_SEQ.resize(BATCH_SIZE+2);
	PE2_QUAL.resize(BATCH_SIZE+2);

	int len_pe1;
	int len_pe2;

	gzFile fp_pe1;
	kseq_t *seq_pe1;

	gzFile fp_pe2;
	kseq_t *seq_pe2;

	fp_pe1 = gzopen((P2In->InFq1).c_str(), "r");
	seq_pe1 = kseq_init(fp_pe1);

	fp_pe2 = gzopen((P2In->InFq2).c_str(), "r");
	seq_pe2 = kseq_init(fp_pe2);

	string id_1, seq_1, plus_1, qual_1 ;
	string id_2, seq_2, plus_2, qual_2 ;

	string OUT=(P2In->OutPrefix);

	string outname_pe1=OUT+"_pe_1_tmp.fa";
	string outname_pe2=OUT+"_pe_2_tmp.fa";

	if ((P2In->MinCount)==1){
		outname_pe1=OUT+"_pe_1.fa";
		outname_pe2=OUT+"_pe_2.fa";
	}

	if (!P2In->OutFa) {
		if ((P2In->MinCount)==1){
			outname_pe1=OUT+"_pe_1.fq";
			outname_pe2=OUT+"_pe_2.fq";
		}else{
			outname_pe1=OUT+"_pe_1_tmp.fq";
			outname_pe2=OUT+"_pe_2_tmp.fq";
		}
	}

	FilePath.push_back(outname_pe1);
	FilePath.push_back(outname_pe2);

	int hq_number=0;

	ofstream OUTH_PE1,OUTH_PE2;
	OUTH_PE1.open(outname_pe1.c_str());
	OUTH_PE2.open(outname_pe2.c_str());
	int seq_num=0;
	int A=0;
	for (A=0 ; A<(P2In->ReadNumber) && ((len_pe1 = kseq_read(seq_pe1)) >= 0) && 
										((len_pe2 = kseq_read(seq_pe2)) >= 0); A++){
		PE1_ID[seq_num]=seq_pe1->name.s;
		PE1_SEQ[seq_num]=seq_pe1->seq.s;
		PE1_QUAL[seq_num]=seq_pe1->qual.s;
		
		PE2_ID[seq_num]=seq_pe2->name.s;
		PE2_SEQ[seq_num]=seq_pe2->seq.s;
		PE2_QUAL[seq_num]=seq_pe2->qual.s;

		seq_num++;
		if (seq_num==BATCH_SIZE) {
			int out_number=Run_PE_low_qual_batch(P2In, seq_num, OUTH_PE1, OUTH_PE2, 
									PE1_ID, PE1_SEQ, PE1_QUAL,
									PE2_ID, PE2_SEQ, PE2_QUAL);
			hq_number+=out_number;
			
			seq_num=0;
		}
	}

	//filter the remaining sequence and output
	if (seq_num!=0) {
		int out_number=Run_PE_low_qual_batch(P2In, seq_num, OUTH_PE1, OUTH_PE2, 
								PE1_ID, PE1_SEQ, PE1_QUAL,
								PE2_ID, PE2_SEQ, PE2_QUAL);
		hq_number+=out_number;
		seq_num=0;
	}

	OUTH_PE1.close();
	OUTH_PE2.close();

	if (A==(P2In->ReadNumber)){
		if  (!((len_pe1 = kseq_read(seq_pe1)) >= 0)) {
			cout<<"INFO: ALL reads "<<A*2<<" are read done"<<endl;
		}
	} else {
		cout<<"INFO: ALL reads "<<A*2<<" are read done"<<endl;
	}

	if ((P2In->MinCount)==1){
		cout<<"INFO: output PE1 read number is "<<hq_number<<"\n";
		cout<<"INFO: output PE2 read number is "<<hq_number<<"\n";
		cout<<"INFO: paired PE1 read number is "<<hq_number<<"\n";
		cout<<"INFO: paired PE2 read number is "<<hq_number<<"\n";
	}

	kseq_destroy(seq_pe1);
	gzclose(fp_pe1);

	kseq_destroy(seq_pe2);
	gzclose(fp_pe2);

	return 0;
}

int Run_SE_low_qual_filter (Para_A24 * P2In, vector<std::string> & FilePath) {

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  SE_ID ;
	vector <string>  SE_SEQ ;
	vector <string>  SE_QUAL ;
	
	SE_ID.resize(BATCH_SIZE+2);
	SE_SEQ.resize(BATCH_SIZE+2);
	SE_QUAL.resize(BATCH_SIZE+2);

	int len_se;
	gzFile fp_se;
	kseq_t *seq_se;

	fp_se = gzopen((P2In->InFq1).c_str(), "r");
	seq_se = kseq_init(fp_se);
	
	string id, seq, plus, qual;

	string OUT=(P2In->OutPrefix);
	string outname_se=OUT+"_tmp.fa";

	if ((P2In->MinCount)==1){
		outname_se=OUT+".fa";
	}

	if (!P2In->OutFa) {
		if ((P2In->MinCount)==1){
			outname_se=OUT+".fq";
		}else{
			outname_se=OUT+"_tmp.fq";
		}
	}

	FilePath.push_back(outname_se);

	ofstream OUTH_SE;
	OUTH_SE.open(outname_se.c_str());

	int hq_number=0;
	int seq_num=0;
	int A=0;

	for (A=0 ; A<(P2In->ReadNumber) && ((len_se = kseq_read(seq_se)) >= 0); A++){
		SE_ID[seq_num]=seq_se->name.s;
		SE_SEQ[seq_num]=seq_se->seq.s;
		SE_QUAL[seq_num]=seq_se->qual.s;
		
		seq_num++;
		if (seq_num==BATCH_SIZE) {
			int out_number=Run_SE_low_qual_batch(P2In, seq_num, OUTH_SE, SE_ID, SE_SEQ, SE_QUAL);
			hq_number+=out_number;
			seq_num=0;
		}
	}

	if (seq_num!=0) {
		int out_number=Run_SE_low_qual_batch(P2In, seq_num, OUTH_SE, SE_ID, SE_SEQ, SE_QUAL);
		hq_number+=out_number;
		seq_num=0;
	}

	OUTH_SE.close();

	if (A==(P2In->ReadNumber)){
		if (!((len_se = kseq_read(seq_se)) >= 0)) {
			cout<<"INFO: ALL reads "<<A<<" are read done"<<endl;
		}
	} else {
		cout<<"INFO: ALL reads "<<A<<" are read done"<<endl;
	}

	if ((P2In->MinCount)==1){
		cout<<"INFO: output SE read number is "<<hq_number<<endl;
	}

	kseq_destroy(seq_se);
	gzclose(fp_se);

	return 0;
}
int GetMinReadKmerCount(Para_A24 * P2In, int read_length){
	int ReadKmerCount=int((read_length-(P2In->Kmer))/(P2In->Windows))+1;
	int MinReadKmerCount=int(ReadKmerCount * (P2In->MinReadKmerRatio));
	return MinReadKmerCount;
}
void GetMinCount(Para_A24 * P2In, const kc_c4x_t *h){
	if(((P2In->MinCount)==0) || (P2In->KmerStatOut)){
		hist_aux_t a;
		uint64_t cnt[256];
		int i, j;
		a.h = h;

		CALLOC(a.cnt, n_thread);
		kt_for(n_thread, worker_hist, &a, 1<<h->p);
		for (i = 0; i < 256; ++i) { cnt[i] = 0;}
		for (j = 0; j < n_thread; ++j) {
			for (i = 0; i < 256; ++i) {
				cnt[i] += a.cnt[j].c[i];
			}
		}
		free(a.cnt);

		uint64_t  Sum=0;
		for (i = 1; i < 256; ++i) {
			Sum+=cnt[i]*i;
		}
		uint64_t  Half=Sum/2;

		if ((P2In->KmerStatOut)) {
			string Kmer_out=(P2In->OutPrefix)+".KmerFre.txt";
			ofstream OUTH_kmer;
			OUTH_kmer.open(Kmer_out.c_str());
			for (i = 1; i < 256; ++i) {
				OUTH_kmer<<i<<"\t"<<cnt[i]<<"\n";
			}
			OUTH_kmer.close();
		}
		else {
			Sum=0;
			for (i = 1; i < 256 &&  Sum<Half; ++i) {
				Sum+=cnt[i]*i;
			}

			if ((P2In->MinCount)==0) { (P2In->MinCount)=i;}
		}
	}

	if ((P2In->MinCount)<2) {P2In->MinCount=2;}
	cout<<"INFO: min kmer count is set to :"<<(P2In->MinCount)<<endl;
}

int Filter_PE_low_kmer_reads(Para_A24 * P2In, const kc_c4x_t *h, bool * PASS_PE1, bool * PASS_PE2, 
			int & Start, int & End, vector <string> & PE1_SEQ, vector <string> & PE2_SEQ){

	int Count_PE1=0; // count of high freq kmer in one read
	int Count_PE2=0;
	
	for (int ii=Start; ii<End; ii++) {
		Count_PE1=ReadHitNum(h, P2In->Kmer, P2In->Windows, P2In->MinCount, PE1_SEQ[ii]);
		Count_PE2=ReadHitNum(h, P2In->Kmer, P2In->Windows, P2In->MinCount, PE2_SEQ[ii]);
		int pe1_length=PE1_SEQ[ii].length();
		int pe2_length=PE2_SEQ[ii].length();
		int pe1_min_count=GetMinReadKmerCount(P2In,pe1_length);
		int pe2_min_count=GetMinReadKmerCount(P2In,pe2_length);
		if (Count_PE1 < pe1_min_count) {
			PASS_PE1[ii]=false;
		}
		else {
			PASS_PE1[ii]= true;
		}

		if (Count_PE2 < pe2_min_count) {

			PASS_PE2[ii]=false;
		}
		else {
			PASS_PE2[ii]= true;
		}

	}

	return 1;
}

void Filter_SE_low_kmer_reads(Para_A24 * P2In, const kc_c4x_t *h, bool * PASS, int &Start, int &End, vector <string> &SEQ) {
	int Count=0;
	for (int ii=Start; ii<End; ii++){
		Count=ReadHitNum(h, P2In->Kmer, P2In->Windows, P2In->MinCount, SEQ[ii]);
		int length=SEQ[ii].length();
		int min_count=GetMinReadKmerCount(P2In, length);
		if (Count < min_count){

			PASS[ii]=false;
		}
		else {
			PASS[ii]= true;
		}
	}
}

vector<int> Run_PE_low_kmer_batch(Para_A24 * P2In, const kc_c4x_t *h, int &seq_num, 
				ofstream &OUTH_PE1, ofstream &OUTH_PE2, ofstream &OUTH_SE1, ofstream &OUTH_SE2,
				vector <string> &PE1_ID, vector <string> &PE1_SEQ, vector <string> &PE1_QUAL,
				vector <string> &PE2_ID, vector <string> &PE2_SEQ, vector <string> &PE2_QUAL){
	
	std::vector<std::thread> threads;
	int * Start =new int [n_thread];
	int * End =new int [n_thread];
	bool  *PASS_PE1 =new bool [BATCH_SIZE];
	bool  *PASS_PE2 =new bool [BATCH_SIZE];

	for (int i = 0; i < n_thread; i++){
		Start[i]=i*BinWind;
		End[i]=Start[i]+BinWind;
		if (End[i]>seq_num) {End[i]=seq_num;} if (Start[i]>=End[i]) {continue;}
		threads.push_back(thread(Filter_PE_low_kmer_reads,P2In,h,PASS_PE1,PASS_PE2,
								 ref(Start[i]),ref(End[i]),
						         ref(PE1_SEQ),ref(PE2_SEQ)));
	}

	for (auto& thread : threads) {
		thread.join();
	}

	threads.clear();
	
	vector<int> hfk_number=Out_PE_seq(P2In, seq_num, OUTH_PE1, OUTH_PE2, OUTH_SE1, OUTH_SE2, 
			   PASS_PE1, PASS_PE2, PE1_ID, PE1_SEQ, PE1_QUAL, PE2_ID, PE2_SEQ, PE2_QUAL);

	delete [] Start;
	delete [] End;
	delete [] PASS_PE1;
	delete [] PASS_PE2;

	return hfk_number;
}

int Run_SE_low_kmer_batch(Para_A24 * P2In, const kc_c4x_t *h, int &seq_num, ofstream &OUTH_SE, 
				vector <string> &SE_ID, vector <string> &SE_SEQ, vector <string> &SE_QUAL){
	
	std::vector<std::thread> threads;
	int * Start =new int [n_thread];
	int * End =new int [n_thread];
	bool  *PASS =new bool [BATCH_SIZE];

	for (int i = 0; i < n_thread; i++){
		Start[i]=i*BinWind;
		End[i]=Start[i]+BinWind;
		if (End[i]>seq_num) {End[i]=seq_num;} if (Start[i]>=End[i]) {continue;}
		threads.push_back(thread(Filter_SE_low_kmer_reads,P2In,h,PASS,
								 ref(Start[i]),ref(End[i]), ref(SE_SEQ)));
	}

	for (auto& thread : threads) {
		thread.join();
	}

	threads.clear();

	int hfk_number=Out_SE_seq(P2In, seq_num, OUTH_SE, PASS, SE_ID, SE_SEQ, SE_QUAL);

	delete [] Start;
	delete [] End;
	delete [] PASS;

	return hfk_number;
}

int Run_PE_low_kmer_filter(Para_A24 * P2In, vector<std::string> &FilePath, const kc_c4x_t *h){
	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  PE1_ID;
	vector <string>  PE1_SEQ;
	vector <string>  PE1_QUAL;
	
	vector <string>  PE2_ID ;
	vector <string>  PE2_SEQ ;
	vector <string>  PE2_QUAL ;

	PE1_ID.resize(BATCH_SIZE+2);
	PE1_SEQ.resize(BATCH_SIZE+2);
	PE1_QUAL.resize(BATCH_SIZE+2);
	
	PE2_ID.resize(BATCH_SIZE+2);
	PE2_SEQ.resize(BATCH_SIZE+2);
	PE2_QUAL.resize(BATCH_SIZE+2);
	
	int len_pe1;
	int len_pe2;

	gzFile fp_pe1;
	kseq_t *seq_pe1;

	gzFile fp_pe2;
	kseq_t *seq_pe2;

	fp_pe1 = gzopen(FilePath[0].c_str(), "r");
	seq_pe1 = kseq_init(fp_pe1);

	fp_pe2 = gzopen(FilePath[1].c_str(), "r");
	seq_pe2 = kseq_init(fp_pe2);

	string id_1, seq_1, plus_1, qual_1 ;
	string id_2, seq_2, plus_2, qual_2 ;

	string OUT=(P2In->OutPrefix);
	string outname_pe1=OUT+"_pe_1.fa";
	string outname_pe2=OUT+"_pe_2.fa";
	string outname_se1=OUT+"_se_1.fa";
	string outname_se2=OUT+"_se_2.fa";

	if (!P2In->OutFa) {
		outname_pe1=OUT+"_pe_1.fq";
		outname_pe2=OUT+"_pe_2.fq";
		outname_se1=OUT+"_se_1.fq";
		outname_se2=OUT+"_se_2.fq";
	}

	int hfk_pe_number=0;
	int hfk_se_number=0;
	ofstream OUTH_PE1,OUTH_PE2,OUTH_SE1,OUTH_SE2;
	OUTH_PE1.open(outname_pe1.c_str());
	OUTH_PE2.open(outname_pe2.c_str());
	OUTH_SE1.open(outname_se1.c_str());
	OUTH_SE2.open(outname_se2.c_str());
	int seq_num=0;
	int A=0;
	
	for (A=0 ; A<(P2In->ReadNumber) && ((len_pe1 = kseq_read(seq_pe1)) >= 0) && 
										((len_pe2 = kseq_read(seq_pe2)) >= 0); A++){
		PE1_ID[seq_num]=seq_pe1->name.s;
		PE1_SEQ[seq_num]=seq_pe1->seq.s;
		PE1_QUAL[seq_num]=seq_pe1->qual.s;
		
		PE2_ID[seq_num]=seq_pe2->name.s;
		PE2_SEQ[seq_num]=seq_pe2->seq.s;
		PE2_QUAL[seq_num]=seq_pe2->qual.s;

		seq_num++;
		if (seq_num==BATCH_SIZE) {
			vector<int> out_number=Run_PE_low_kmer_batch(P2In, h, seq_num, 
									OUTH_PE1, OUTH_PE2, OUTH_SE1, OUTH_SE2, 
									PE1_ID, PE1_SEQ, PE1_QUAL,
									PE2_ID, PE2_SEQ, PE2_QUAL);
			hfk_pe_number+=out_number[0];
			hfk_se_number+=out_number[1];
			seq_num=0;
		}
	}

	//filter the remaining sequence and output
	if (seq_num!=0) {
		vector<int> out_number=Run_PE_low_kmer_batch(P2In, h, seq_num, 
								OUTH_PE1, OUTH_PE2, OUTH_SE1, OUTH_SE2, 
								PE1_ID, PE1_SEQ, PE1_QUAL,
								PE2_ID, PE2_SEQ, PE2_QUAL);
		hfk_pe_number+=out_number[0];
		hfk_se_number+=out_number[1];
		seq_num=0;
	}

	OUTH_PE1.close();
	OUTH_PE2.close();
	OUTH_SE1.close();
	OUTH_SE2.close();

	cout<<"INFO: output PE1 read number is "<<hfk_pe_number+hfk_se_number<<"\n";
	cout<<"INFO: output PE2 read number is "<<hfk_pe_number+hfk_se_number<<"\n";
	cout<<"INFO: paired PE1 read number is "<<hfk_pe_number<<"\n";
	cout<<"INFO: paired PE2 read number is "<<hfk_pe_number<<"\n";
	cout<<"INFO: un-paired PE1 read number is "<<hfk_se_number<<"\n";
	cout<<"INFO: un-paired PE2 read number is "<<hfk_se_number<<"\n";

	kseq_destroy(seq_pe1);
	gzclose(fp_pe1);

	kseq_destroy(seq_pe2);
	gzclose(fp_pe2);

	return 0;
}

int Run_SE_low_kmer_filter( Para_A24 * P2In, vector<std::string> & FilePath, const kc_c4x_t * h){
	
	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  SE_ID ;
	vector <string>  SE_SEQ ;
	vector <string>  SE_QUAL ;
	
	SE_ID.resize(BATCH_SIZE+2);
	SE_SEQ.resize(BATCH_SIZE+2);
	SE_QUAL.resize(BATCH_SIZE+2);

	int len_se;
	gzFile fp_se;
	kseq_t *seq_se;

	fp_se = gzopen(FilePath[0].c_str(), "r");
	seq_se = kseq_init(fp_se);
	
	string id, seq, plus, qual;

	string OUT=(P2In->OutPrefix);
	string outname_se=OUT+".fa";

	if (!P2In->OutFa) {
		outname_se=OUT+".fq";
	}

	int hfk_number=0;

	ofstream OUTH_SE;
	OUTH_SE.open(outname_se.c_str());

	int seq_num=0;
	int A=0;

	for (A=0 ; A<(P2In->ReadNumber) && ((len_se = kseq_read(seq_se)) >= 0); A++){
		SE_ID[seq_num]=seq_se->name.s;
		SE_SEQ[seq_num]=seq_se->seq.s;
		SE_QUAL[seq_num]=seq_se->qual.s;
		
		seq_num++;
		if (seq_num==BATCH_SIZE) {
			int out_number=Run_SE_low_kmer_batch(P2In, h, seq_num, OUTH_SE, SE_ID, SE_SEQ, SE_QUAL);
			hfk_number+=out_number;
			seq_num=0;
		}
	}

	if (seq_num!=0) {
		int out_number=Run_SE_low_kmer_batch(P2In, h, seq_num, OUTH_SE, SE_ID, SE_SEQ, SE_QUAL);
		hfk_number+=out_number;
		seq_num=0;
	}

	OUTH_SE.close();

	cout<<"INFO: output SE read number is "<<hfk_number<<endl;

	kseq_destroy(seq_se);
	gzclose(fp_se);

	return 0;
}

string GetFileExtension(const string& FilePath) {
	size_t dotPos = FilePath.rfind('.');
	if (dotPos == string::npos) {
		return "";
	} 
	else {
		return FilePath.substr(dotPos + 1);
	}
}

void Check_outprefix(Para_A24 * P2In) {
	string path=P2In->OutPrefix;
	string ext = GetFileExtension(path);

	if (ext == "gz") {
		path = path.substr(0, path.rfind('.'));
		P2In->OutPrefix=path;
		ext = GetFileExtension(path);
	}
	if((ext == "fq") || (ext == "fastq") || (ext == "fa") || (ext == "fasta")){
		path = path.substr(0, path.rfind('.'));
		P2In->OutPrefix=path;
	}
}

//////////////////main///////////////////
int main (int argc, char *argv[]) {
	Para_A24 * P2In = new Para_A24;
	int InPESE=1; // 1 for PE; 2 for SE; 0 for Unknow
	InPESE=parse_cmd_HFKreads(argc, argv, P2In);
	if(InPESE==0) {
		delete P2In ;
		return 0 ;
	}

	Check_outprefix(P2In);

	if  (InPESE==2) {
		(P2In->InFq1)=P2In->InSeFq; // use for qtype 
	}

	P2In->MaxQV=Get_qType((P2In->InFq1),P2In); // 0 for fasta

	int read_length=(P2In->ReadLength);
	cout<<"INFO: middle read length is "<<read_length<<" bp"<<endl;
	
	if (VECMAX == 102400){
		if (InPESE==1){
			VECMAX=(VECMAX*100)/(read_length*2);
		} else if (InPESE==2){
			VECMAX=(VECMAX*100)/read_length;
		}
	}
	
	BinWind = ceil(VECMAX/n_thread);
	BATCH_SIZE = BinWind*n_thread;

	for (int i=0; i<256; i++) {
		NArry[i]=false;
		LowArry[i]=true;
	}
	NArry['N']=true ; NArry['n']=true ;

	for (int i=0; i<(P2In->MinBaseQ); i++) {
		LowArry[i]=false;
	}

	vector<std::string> FilePath;
	if (InPESE==1){
		P2In->ReadNumber=(P2In->ReadNumber)/2;
		Run_PE_low_qual_filter(P2In,FilePath);
	} else if(InPESE==2){
		Run_SE_low_qual_filter(P2In,FilePath);
	}

	if ((P2In->MinCount)==1) {
		delete P2In ;
		return (0);
	}

	//get kmer freq 
	uint64_t hash_size = 10000000; // Initial size of hash.
	if (n_thread>10) {
		hash_size=int(n_thread*hash_size/10);
	}

	kc_c4x_t *h;
	int p=10;  // Minimum length of counting field

	// create the hash
	h = count_file(FilePath, P2In->Kmer, P2In->Windows, p, hash_size, n_thread);

	GetMinCount(P2In,h);
	
	if (InPESE==1){
		Run_PE_low_kmer_filter(P2In, FilePath, h);
	} else if(InPESE==2){
		Run_SE_low_kmer_filter(P2In, FilePath, h);
	}

	for (auto& file : FilePath){
		remove(file.c_str());
	}

	for (int i = 0; i < 1<<p; ++i) {
		kc_c4_destroy(h->h[i]);
	}

	free(h->h); free(h);

	delete P2In ;
	return (0);
}

#endif  //
