# Iterative-blast-v.1.0
The codes and tools for iterative blast (iblast)

Please first download blast program from NCBI and install, and then put the .exe or .py tool files into the bin directory of the installed blast program before running on Windows system or in Python environment, respectively.

Four blast tools are included: UTR-iblastn, iblastx, iblastp, and iblastn/itblastx.

In the input box for blast parameters, please note that:
1)	format of the intermediate output file is fixed as -outfmt 6, please ignore this parameter;
2)	type, for example, -evalue 0.0001, for parameter of threshold e-value, according to the NCBI blast;
3)	other available parameter items and their values are -word_size <Integer, >=4>, -gapopen <Integer>, -gapextend <Integer>, -penalty <Integer, <=0>, -reward <Integer, >=0>, -num_threads <Integer, >=1>, check their usages in the NCBI blast program, and use a parameter item by typing in like -num_threads 2; 
4)	other blast parameter items are not included, make your own extensions if necessary.

Test files are attached in the Test files.zip, which are used for tool running test.

The final results from tool running include:
1)	B_total.text, which is a collection of annotated sequences from query file against db file in the iterative blast analysis;
2)	A_iterations.txt, indicating number of times (make this n) that a blast analysis iterated;
4)	<query_file_name> n.fa, the remaining sequences not annotated by iterative blast with the n iteration.

Please cite "Zhang S et al. 2023. Conserved untranslated regions of multipartite viruses: Natural markers of novel viral genomic components and tags of viral evolution"
