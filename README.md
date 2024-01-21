# PanDepth
<b>HFKR: An ultra-fast and efficient tool for high frequency kmer sequencing reads extraction</b>

##  1. Install
### (1) Pre-built binaries for x86_64 Linux
```
wget -c https://github.com/HuiyangYu/HFKR/releases/download/v2.00/HFKR-2.00-Linux-x86_64.tar.gz
tar zvxf HFKR-2.00-Linux-x86_64.tar.gz
cd HFKR-2.00-Linux-x86_64
./hfkreads -h
```
### (2) Building from source
```
git clone https://github.com/HuiyangYu/HFKR.git
cd HFKR
make
cd bin
./hfkreads -h
```
## 2. Usage
```
Usage: hfkreads -1 PE1.fq.gz -2 PE2.fq.gz -o OutFrefix
 Input/Output options:
   -1	<str>   paired-end fasta/q file1
   -2	<str>   paired-end fasta/q file2
   -s	<str>   single-end fasta/q
   -o	<str>   prefix of output file
 Filter options:
   -b	<int>   min base quality [0]
   -q	<int>   min average base quality [20]
   -l	<int>   min length of read [half]
   -r	<float> max unknown base (N) ratio [0.1]
   -k	<int>   kmer length [31]
   -w	<int>   window size [10]
   -m	<int>   min kmer count for high freq kmer [2]
   -x	<int>   min count of read with high freq kmer [10]
   -n	<int>   read number to use [1000000]
   -a	        use all the read number
 Other options:
   -c           compress the outPut File
   -f           outPut the KmerFre File
   -A           keep output quality info
   -D           force PE SE pair-wise File 
   -t           thread to run [4]
   -h           show help [v2.00]
```
## 3. Example
### 3.1 Extracting high-frequency k-mer reads from paired-end sequencing reads
```
hfkreads -1 PE_1.fq.gz -2 PE_2.fq.gz -o test1
```
The output files consist of four files: test1_pe_1.fa, test1_pe_2.fa, test1_se_1.fa, and test1_se_2.fa. The reads with the 'se' label are unpaired high-frequency k-mer reads.

### 3.2 Extracting high-frequency k-mer reads from single-end sequencing reads
```
hfkreads -s SE.fq.gz -o test2
```
The output result consists of a single file named "test2.fa".

### 3.3 Extracting high-quality reads without filtering the k-mers
```
hfkreads -1 PE_1.fq.gz -2 PE_2.fq.gz -m 1 -o test3
hfkreads -s SE.fq.gz -m 1 -o test4
```
To filter out low-quality reads only, the default '-m 2' parameter should change to '-m 1'.

### 3.4 Other options
The '-A' parameter outputs quality values in fastq format, while the default is to output files in fasta format. 
The '-c' parameter compresses the output file, and by default, the output file is uncompressed.

## 4. License
-------

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
