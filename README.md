# Iterative-blast-v.1.0
The codes and tools for iterative blast (iblast)

The environments, Python module, and independent open-source program for tests:
Windows, Version=11;
Python, Version=3.10.2;
ttkbootstrap in Python, Version=1.10.1;
BLAST, Version=2.12.0+.

Please first download blast program from NCBI and install, and then put the .exe or .py tool files into the bin directory of the installed blast program before running on Windows system or in Python environment, respectively.

Four blast tools are included: UTR-iblastn, iblastx, iblastp, and iblastn/itblastx.

In the input box for blast parameters, please note that:
1)	format of the intermediate output file is fixed as -outfmt 6, please ignore this parameter;
2)	type, for example, -evalue 0.0001, for parameter of threshold e-value, according to the NCBI blast;
3)	other available parameter items and their values, please check their usages in the NCBI blast program.

Test files are attached for each tool, which are used for running test.

The final results from tool running include:
1)	B_total.text, which is a collection of annotated sequences from query file against db file in the iterative blast analysis;
2)	A_iterations.txt, indicating number of times (make this n) that a blast analysis iterated;
4)	<query_file_name>n.fa, the remaining sequences not annotated by iterative blast with the n iteration.
