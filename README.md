# Iterative-blast-v.1.0
The program and codes iterative blast (iblast)

Please first download blast program from NCBI and install, and then put the .exe or .py tool files into the bin directory of the installed blast program before running on Windows system or in Python environment, respectively.

Four blast tools are included: UTR-iblastn, iblastx, iblastp, and iblastn/itblastx.

In the input box for blast parameters, please note that:
1)	format of the intermediate output file is fixed as -outfmt 6, and this file is unavailable during and after the analysis;
2)	for -evalue parameter item, just manually type in a value being <Real> (>0), such as e-4 instead of -evalue e-4;
3)	other available parameter items and their values are -word_size <Integer, >=4>, -gapopen <Integer>, -gapextend <Integer>, -penalty <Integer, <=0>, -reward <Integer, >=0>, -num_threads <Integer, >=1>, check their usages in the NCBI blast program, and use a parameter item by typing in like -num_threads 2; 
4)	other blast parameter items are not included, make your own extensions if necessary.

Test files are attached in the Test files.zip, which are used for tool running test.

Default input file names for query and database (db) building are A.fasta.fa and B.fasta.fa respectively, and the final results from tool running include:
1)	B_total.text, which is a collection of annotated A.fasta.fa sequences from iterative blast analysis;
2)	A_iterations.fasta.fa.txt, which is iteration time (make this n) of blast analysis;
3)	An.fasta.fa, the remaining sequences not annotated in A.fasta.fa by iterative blast analysis.

Please cite "Zhang S et al. 2023. Conserved untranslated regions of multipartite viruses: Natural markers of novel viral genomic components and tags of viral evolution. Virus Evolution."
