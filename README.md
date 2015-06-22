MAPGD version 2.2

<h5> To download this program click the "Download ZIP" button to the upper right -> </h5>

Copyright (C) Michael Lynch, see notice at end of README. 

<h3> Introduction </h3>

<h5> Important changes from 2.0 </h5>

Changes to the interface have been made to better conforming to POSIX guidlines:
The the syntax of the -p option (which specifies the populations to be analized) has been changed from a space delimited list to a comma delimited list. See the -h option from more information.
The default behavior of the cp and ep commands bas been changed to read and write to the stdin and stdout (respectively).

-i and -o options will continue to function as they have previously. 

Headers have been added to .pro files to preserve information about input files and aid with compatibility as development continues. All .pro files w/o headers are assumed to have been generated by v 2.0 with the -c 6 option.

<h5> the ep and cp commands </h5>

These commands use a maximum-likelihood (ML) procedure to estimate allele-frequencies from the numbers of the four nucleotides (quartets) observed at individual genomic sites. For each site, the major and minor nucleotides are identified, their frequencies are estimated by ML.

The ep command estimates whether a site is significantly polymorphic within a population; The cp command estimates whether the frequency of major alleles at a site differ between two sample population.

The designated major nucleotide is simply the one with the highest rank, and the minor nucleotide is the one with the second highest rank. If the top three ranks are all equal, the site is treated as unresolvable, with both major and minor nucleotides designated in the output by a *.

If the second and third ranks are equal but lower than the major-nucleotide count, the site is treated as monomorphic, with the minor nucleotide again designated by a *.

The ep command is useful to distinguish between sites that may appear polymorphic because of sequencing errors and sites that are truly polymorphic, and the cp command is useful to determine whether sampeling and sequencing errors can explain differences in estimates of allele frequencies, or whether allele frequencies differ between population samples.
<h5> options : </h5>

<h3> Making the Input File </h3>

<h5> the .pro File </h5>
A number of .pro File formats are excepted by mapgd.

The input file is a plain text file consists of five tab delimited columns, one for each site: the first entry is an arbitrary identifier (e.g., and site), and the final four are integer values for the number of times an A, C, G, and T was observed at the site. 

We call this file format the .pro file format. Files in this format can be generated from mpileup files using the command "mapgd proview" 

Currently mapgd allows for the estimation of allele frequencies in pooled population sequence through the "ep" command, the comparison of allele frequencies between two populations with the "cp" command, and the conversion of mpileups to the ".pro" format though the proview command. We hope to add more commands in the future.

For example, if the sequencing center gives you two files called "population1.fastq" and "population1.fastq" your entire work flow might look something like this: 

<h5> a typical workflow </h5>
	bwa aln Reference.fna population1.fastq > population1.sai
	bwa sam Reference.fna population1.sai population1.fastq > reads.sam
	samtools view -bS population1.sam > population1.bam
	samtools sort population1.bam population1.sort
	samtools index population1.sort.bam

	bwa aln Reference.fna population2.fastq > population2.sai
	bwa sam Reference.fna population2.sai population2.fastq > population2.sam
	samtools view -bS population2.sam > population2.bam
	samtools sort population2.bam population2.sort
	samtools index population2.sort.bam

	samtools mpileup -q 25 -Q 25 -B population1.sort.bam population2.sort.bam > metapopulation.mpileup
	mapgd proview -i metapopulation.pileup > metapopulation.pro

	mapgd ep -i metapopulation.pro -p 1,2 -o allelfrequency_estimates.txt
	mapgd cp -i metapopulation.pro -p 1,2 -o allelfrequency_comparison.txt


Alternatevely, you may not desire to create .mpileup or .pro files, in which case you can perform an analysis more quickly using the following syntax:

<h5> the work flow using  I/O redirection </h5>
        bwa aln Reference.fna population1.fastq > population1.sai
        bwa sam Reference.fna population1.sai population1.fastq > reads.sam
        samtools view -bS population1.sam > population1.bam
        samtools sort population1.bam population1.sort
        samtools index population1.sort.bam

        bwa aln Reference.fna population2.fastq > population2.sai
        bwa sam Reference.fna population2.sai population2.fastq > population2.sam
        samtools view -bS population2.sam > population2.bam
        samtools sort population2.bam population2.sort
        samtools index population2.sort.bam

        samtools mpileup -q 25 -Q 25 -B population1.sort.bam population2.sort.bam > metapopulation.mpileup | mapgd proview | mapgd ep -p 1,2 > allelefrequency_estimates.txt
        samtools mpileup -q 25 -Q 25 -B population1.sort.bam population2.sort.bam > metapopulation.mpileup | mapgd proview | mapgd cp -p 1,2 > allelefrequency_comparison.txt

<h5> Additional Programs </h5>

This program is intended for use with ".pro" described in the previous section. A slightly modified version of the program sam2pro, written by Bernhard Haubold, is included in this package and can be run by typing "mapgd proview". Programs for converting other file formats to ".pro" files are available in stand alone form http://guanine.evolbio.mpg.de/mlRho/

The workflow described above requires the programs bwa and samtools.

To download bwa please visit http://bio-bwa.sourceforge.net/

To download samtools please visit http://www.htslib.org/

<h5> For windows users </h5>

The windows binaries for V 2.0 are currently unavalible. You can try to compile the source code yourself, as I have tried to refrain from using any platform specific libraries, or you can send me an e-mail telling me to get the binaries up ASAP. 

After clicking the "Download ZIP" button you will be prompted to save or open the file MAPGD-master.zip. Extract this file to the directory of your choice, which should create the new directory "MAPGD-master". Open this directory in windows explorer and click on the file "RunMeWin-32.bat" if you have a 32 bit version of windows or the "RunMeWind-64.bat" if you have a 64 bit verions of windows. If you don't know what version of windows you have, then just try clicking on both. You will be prompted to enter the name of the file you wish to analyze. You can type "test\test.pileup" to analyze test data. To analyze your own data, just drag and drop the file into MAPGD-master folder, click on the RunMeWin-32.bat or RunMeWin-64.bat and type the filename of the file you would like to analyze. Output will be saved to the file output.txt.   

mapgd can also be run from the command prompt, which can be accessed pressing the Windows logo key and r key simultaneously then typing "cmd" into menu which appears.

Although mapgd can be run by windows users, bwa, samtools, and other bioinformatic programs may not be availble in windows.

<h5> Linux or Mac users </h5>


After clicking the "Downlaod ZIP" button you will be prompted to save or open the file MAPGD-master.zip. Save this this file to the directory of your choice, then go to this director in a terminal, for example “cd /home/LynchLab/Downloads/” 

Then type:

	unzip MAPGD.zip
	cd MAPGD-master/src
	make

The program can be installed for all users of a computer by typing:

	sudo make install

Scripts for both Linux and Mac users are present in the top level directory, or the program can be run by typing "mapgd ep -i FILENAME" where FILENAME is the a .pro file.

<h3> The output of ep </h3>

Columns 1 and 2 are site identifiers (ID1 and ID2); 3 and 4 designate major and minor nucleotides (major and minor); 
For each population in the sample four columns are printed : a major allele frequency (freq_P), test statistic for polymorphism (ll_poly), test statistic for fixed for the minor allele (ll_fixed), and coverage (cov);
The final two coloumns are the total coverage (totcov) and the maximum likelihood estimate of the error rate (Error).

Under the assumption of a chi-square distribution for the test statistic with one degree of freedom, significance at the 0.05, 0.01, 0.001 levels requires that the likelihood-ratio test statistic exceed 3.841, 6.635, and 10.827, respectively. 

In principle, the 95% support interval can be obtained by determining the changes in the estimate of the minor allele frequency in both directions required to reduce the log likelihood by the appropriate chi-square value (e.g., 3.841) although this is not currently implemented. 

By default the program prints information to the file "dataout.txt" and this file will appear in the same location as the program. If an alternative file name is desired simply type if "mapgd -o FILENAME" where FILENAME is the name of your output file.

<h3> The output of cp </h3>

Columns 1 and 2 are site identifiers (ID1 and ID2); 3 and 4 designate major and minor nucleotides (major and minor);
For each population two columns are printed : The maximum likelihood estimate of the major allele frequency in that population (freq_P), and the test statistic whether this frequency differes from the average frequency across all samples;

The final two columns are the major allele frequency in the metapopulation (meta_P) and the maximum likelihood estimate of the error rate (Error). Output columns are tab delimited.

By default the program reads information from the stdin and prints information to stdout. File names can be specifide using "mapgd -o FILENAME" where FILENAME is the name of your output file.

<h3> The output of ei </h3>

Columns 1 and 2 are site identifiers (ID1 and ID2); 3 and 4 designate major and minor nucleotides (major and minor); 
A number of other columns exist [TODO : WRITE DESCRIPTION HERE]

Under the assumption of a chi-square distribution for the test statistic with one degree of freedom, significance at the 0.05, 0.01, 0.001 levels requires that the likelihood-ratio test statistic exceed 3.841, 6.635, and 10.827, respectively. 

In principle, the 95% support interval can be obtained by determining the changes in the estimate of the minor allele frequency in both directions required to reduce the log likelihood by the appropriate chi-square value (e.g., 3.841) although this is not currently implemented. 

By default the program prints information to the file "dataout.txt" and this file will appear in the same location as the program. If an alternative file name is desired simply type if "mapgd -o FILENAME" where FILENAME is the name of your output file.

<h3> Reference </h3>

Please cite the following paper when publishing results derived from this program:

Ackerman, M., T. Maruki, and M. Lynch. 2015. MapGD : A program for the maximum likelihood analysis of population-genomic data. ?.

<h3> Copyright Notice </h3>

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

For a copy of the GNU General Public License write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
