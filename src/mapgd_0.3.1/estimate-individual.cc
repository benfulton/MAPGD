/* 

Program estimateIndcpp:

	1) input from a list of site-specific quartets from a labeled population sample;

	2) identify the major and minor alleles, obtain the maximum-likelihood 'like' estimates of allele frequencies, and significance levels.

The designated major nucleotide is simply the one with the highest rank, and the minor nucleotide the one with the second highest rank.
	If the top three ranks are all equal, the site is treated as unresolvable, with both major and minor nucelotides designated in the output by a *.
	If the second and third ranks are equal but lower than the major-nucleotide count, the site is treated as monomorphic, with the minor nucleotide designated by a *.
	
Input File: Three tab delimited dentifiers columns (chromosome, position and ref), followed by an arbitrary number of tab delimited 'quartets', that is counts of the number of times a nucleotide has been observed at a possition.
	Columns are tab delimited, quartets are '/' delimited.
	Default input file is  "datain.txt".

Output File: two columns of site identifiers; reference allele; major allele; minor allele; major-allele frequency; minor-allele frequency; error rate; a ton of other stuff... We really need to clean up the output. 
	Columns are tab delimited.
	Default name is "dataout.txt".
*/

#include "estimate-individual.h"
#include <ciso646>
#include <mpi.h>

#define BUFFER_SIZE 500
#define PRAGMA

enum MESSAGE_TAGS { MAPGD_LINES=1, MAPGD_GOFS };

/*
float_t compare (allele_stat mle1, allele_stat mle2, Locus &site1, Locus &site2,  models &model){
	Locus site3=site1+site2;
	alele_stat mle3;
	maximize_grid(site3, mle3, model, gofs, MIN, MAXGOF, MAXPITCH+texc);
	return mle1.ll+mle2.ll-mle3.ll
}*/

/// Estimates a number of summary statistics from short read sequences.
/**
 */
allele_stat estimate (Locus &site, models &model, std::vector<float_t> &gofs, const count_t &MIN, const float_t &EMLMIN, const float_t &MINGOF, const size_t &MAXPITCH){

	allele_stat mle, temp;					//allele_stat is a basic structure that containes all the summary statistics for
								//an allele. It gets passed around a lot, and a may turn it into a class that has
								//some basic read and write methods.

	mle.gof=0; mle.efc=0; mle.MM=0; mle.Mm=0; mle.mm=0; 	 //Initialize a bunch of summary statics as 0. 
								 //This should be moved over to the constructor of allele_stat 
								 //(when that constructor is writen). I'm a little concerned that
								 //allele_stat has gotten too bloated, but . . . 
	mle.N=0;

	site.mask_low_cov(MIN);
	count_t texc=site.maskedcount(), rexc;
	rexc=texc;

	if (init_params(site, mle, EMLMIN) ){		//If >90% of reads agree, then assume a homozygote,
								//otherwise, assume heterozygote.
	if (mle.null_error!=0){
		rexc=maximize_grid(site, mle, model, gofs, -MINGOF, MAXPITCH+texc);	//trim bad clones and re-fit the model.
//		rexc=maximize_newton(site, mle, model, gofs, MIN, -MINGOF, MAXPITCH+texc);	//trim bad clones and re-fit the model.
	}
	else
		rexc=maximize_analytical(site, mle, model, gofs, -MINGOF, MAXPITCH+texc);	//trim bad clones and re-fit the model.
	}

	mle.excluded=rexc-texc;

	// CALCULATE THE LIKELIHOODS 

	mle.ll=model.loglikelihood(site, mle);		//Sets the site.ll to the log likelihood of the best fit (ll). 

	if (mle.freq<0.5){					//Check to see if the major and minor alleles are reversed.
		std::swap(mle.major, mle.minor);
		std::swap(mle.MM, mle.mm);
		mle.freq=1.-mle.freq;
		site.swap(0, 1);
		
	}
	else if (mle.freq==0.5){
		if (rand() % 2){				//If the major and minor allele frequencies are identical, 
			std::swap(mle.major, mle.minor);	//flip a coin to determine the major and minor allele.
			std::swap(mle.MM, mle.mm);
			mle.freq=1.-mle.freq;
			site.swap(0, 1);
		};
	};

	temp=mle; temp.MM=1.0; temp.Mm=0.; temp.mm=0.;		//Copies site to mono, then sets mono to a monomophic site 
								//(i.e. sets the genotypic frequencies Mm and mm to 0.
	temp.error=mle.null_error;				//Sets the error rate of mono to the null error rate.
	temp.freq=1.;
	temp.f=0.;
	mle.monoll=model.loglikelihood(site, temp);			//Sets the site.ll to the log likelihood of the best fit (ll). 

	if (mle.monoll>mle.ll){
		mle.ll=mle.monoll;
		mle.error=mle.null_error;
	};
	temp=mle; 
	temp.MM=pow(mle.freq, 2);				//Similar set up to mono, but now assuming 
	temp.Mm=2.*mle.freq*(1.-mle.freq); 			//Hardy-Weinberg equilibrium.
	temp.mm=pow(1.-mle.freq, 2);				//?
	mle.hwell=model.loglikelihood(site, temp);		//?
	return mle;
}

void write(std::ostream& out, uint32_t readed, profile& pro, profile& pro_out, allele_stat* buffer_mle, Locus* buffer_site, float MINGOF, float_t A)
{
	for (uint32_t x = 0; x < readed; ++x) {
		// Now print everything to the *out stream, which could be a file or the stdout. 
		//TODO move this over into a formated file.
		//?
		if (2 * (buffer_mle[x].ll - buffer_mle[x].monoll) >= A) {
			out << std::fixed << std::setprecision(6) << pro.getids(buffer_site[x]) << '\t' << buffer_site[x].getname(0) << '\t' << buffer_site[x].getname_gt(1) << '\t';
			out << std::fixed << std::setprecision(6) << buffer_mle[x] << std::endl;
		}
		if (buffer_mle[x].gof < -MINGOF) buffer_site[x].maskall();
		if (pro_out.is_open()) {
			buffer_site[x].id0 = pro_out.encodeid0(pro.decodeid0(buffer_site[x].id0));
			pro_out.write(buffer_site[x]);
		}
	}
}

bool mpi_select(int lineid, int taskid, int num_tasks)
{
	lineid = lineid % (BUFFER_SIZE*num_tasks);
	return lineid >= taskid*BUFFER_SIZE && lineid < (taskid + 1)*BUFFER_SIZE;
}

std::string mpi_recieve_string(int rank)
{
	MPI_Status status;
	MPI_Probe(rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	int count;
	MPI_Get_count(&status, MPI_CHAR, &count);
	char *buf = new char[count];
	MPI_Recv(buf, count, MPI_CHAR, rank, MAPGD_LINES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	std::string result(buf, count);
	delete[] buf;
	return result;
}

void do_estimate(allele_stat* buffer_mle, Locus& buffer_site, models& model, std::vector<int>& ind,
	std::vector <float_t>& sum_gofs, std::vector <float_t>& gofs_read, count_t MIN, float_t EMLMIN, float_t MINGOF,
	count_t MAXPITCH)
{
	std::vector <float_t> gofs(ind.size());
	buffer_site.unmaskall();
	*buffer_mle = estimate(buffer_site, model, gofs, MIN, EMLMIN, MINGOF, MAXPITCH);
	if (2 * (buffer_mle->ll - buffer_mle->monoll) >= 22) {
		for (size_t i = 0; i < sum_gofs.size(); i++) {
			sum_gofs[i] += gofs[i];
			if (gofs[i] != 0) gofs_read[i]++;
		}
	}

}

int estimateInd(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string infile="";
	std::string outfile="";
	std::string outfilepro;

	bool verbose=false;
	bool quite=false;
	bool noheader=false;
	float_t EMLMIN=0.001;
	count_t MIN=0;
	float_t A=0.00;
	float_t MINGOF=2.00;
	count_t MAXPITCH=96;

	count_t skip=0;
	count_t stop=-1;

	std::vector <int> ind;

	/* sets up the help messages and options, see the 'interface.h' for more detials. */

	env_t env;
	env.setname("mapgd ei");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman and Takahiro Maruki");
	env.setdescription("Uses a maximum likelihood approach to estimate population genomic statistics from an individually 'labeled' population.");

	env.optional_arg('i',"input", 	&infile,	&arg_setstr, 	"an error occured while setting the name of the input file.", "the input file for the program (default stdout).");
	env.optional_arg('o',"output", 	&outfile,	&arg_setstr, 	"an error occured while setting the name of the output file.", "the output file for the program (default stdin).");
	env.optional_arg('p',"out-pro", &outfilepro,	&arg_setstr, 	"an error occured while setting the name of the output file.", "name of a 'cleaned' pro file (default none).");
	env.optional_arg('I',"individuals", &ind, 	&arg_setvectorint, "please provide a list of integers", "the individuals to be used in estimates.\n\t\t\t\ta comma seperated list containing no spaces, and the format X-Y can be used to specify a range (defualt ALL).");
	env.optional_arg('m',"minerror", &EMLMIN, 	&arg_setfloat_t, "please provide a float.", "prior estimate of the error rate (defualt 0.001).");

//	env.optional_arg('c',"columns", &EMLMIN, 	&arg_setfloat_t, "please provide a float.", "number of columsn in profile (if applicable).");

	env.optional_arg('M',"mincoverage", &MIN, 	&arg_setint, 	"please provide an int.", "minimum coverage for an individual at a site for an individual to be used (defualt 4).");
	env.optional_arg('a',"alpha", 	&A, 		&arg_setfloat_t, "please provide a float.", "cut-off value for printing polymorphic sites (default 0.0).");
	env.optional_arg('g',"goodfit", &MINGOF,	&arg_setfloat_t, "please provide a float.", "cut-off value for the goodness of fit statistic (defaults 2.0).");
	env.optional_arg('N',"number", 	&MAXPITCH,	&arg_setint, 	"please provide an int.", "cut-off value for number of bad individuals needed before a site is removed entirely (default 96).");
	env.optional_arg('S',"skip", 	&skip,		&arg_setint, 	"please provide an int.", "number of sites to skip before analysis begins (default 0).");
	env.optional_arg('T',"stop", 	&stop,		&arg_setint, 	"please provide an int.", "maximum number of sites to be analyzed (default All sites)");
	env.flag(	'H',"noheader", &noheader,	&flag_set, 	"takes no argument", "disables printing a headerline.");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occured while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occured while displaying the version message.", "prints the program version");
	env.flag(	'V', "verbose", &verbose,	&flag_set, 	"an error occured while enabeling verbose excecution.", "prints more information while the command is running.");
	env.flag(	'q', "quite", 	&quite,		&flag_set, 	"an error occured while enabeling quite execution.", "prints less information while the command is running.");

	if ( parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.

	profile pro, pro_out;		//profile is a fairly complete class that lets us read and write from pro files, 
					//which are files containing set of read 'quartets' that specify the number of 
					//A,C,G and T read at some specific location in a genome. See proFile.h for more info.

	MPI_Init(&argc, &argv);

	int   numtasks, taskid;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	//gcfile out;
	std::ostream *out=&std::cout;
	std::ofstream outFile;

	if (infile.size()!=0) {					//Iff a filename has been set for infile
		pro.open(infile.c_str(), std::fstream::in);	
		if (!pro.is_open() ) {				//try to open a profile of that name.
			printUsage(env);			//Print help message on failure.
		} 
	}
	else {
		pro.open(std::fstream::in);			//Iff no filename has been set for infile, open profile from stdin.
	};

	if (outfile.size()!=0) {
		outFile.open(outfile.c_str(), std::fstream::out);
		if (!outFile.is_open() ) printUsage(env);
		out=&outFile;
	};

	//else out.open('w', CSV);				//Iff no filename has been set for outfile, pgdfile prints to stdout.

	count_t outc=6;
	char cdel='\t';
	char qdel='/';
	bool binary=false;
	if (taskid == 0 && outfilepro.size()!=0) {				//Same sort of stuff for the outf
		if (binary) {
			pro_out.open(outfilepro.c_str(), std::fstream::out | std::fstream::binary);
			if (!pro_out.is_open() ) {printUsage(env); exit(0);}
		} else {
			pro_out.open(outfilepro.c_str(), std::fstream::out);
			if (!pro_out.is_open() ) {printUsage(env); exit(0);}
		}
		pro_out.setsamples(pro.size() );
		pro_out.setcolumns(outc);
		pro_out.set_delim_column(cdel);
		pro_out.set_delim_quartet(qdel);
		for (size_t y=0;y<pro.size(); ++y) pro_out.set_sample_name(y, pro.get_sample_name(y) );
		if (not (noheader) ) pro_out.writeheader();
	};


	/* this is the basic header of our outfile, should probably be moved over to a method in allele_stat.*/
	if (not (noheader) and taskid == 0){
			std::string id1="id1\t";
			switch (pro.get_columns()) {
				case 5:
				case 6:
					*out << id1 << "\tid2\tmajor\tminor\tcov\tM\tm\terror\tnull_e\tf\tMM\tMm\tmm\th\tpoly_ll\thwe_ll\tgof\tef_chrm\tN\tN_cut\tmodel_ll" << std::endl;
				break;
				case 7:
					*out << id1 << "\tid2\tref\tmajor\tminor\tcov\tM\tm\terror\tnull_e\tf\tMM\tMm\tmm\th\tpoly_ll\thwe_ll\tgof\tef_chrm\tN\tN_cut\tmodel_ll" << std::endl;
				break;
			}	
	}
	
	pro.maskall();							//Turn off the ability to read data from all clones by default. 

	if ( ind.size()==0 ) { 						//Iff the vector ind (which should list the clones to 
		ind.clear();						//be read from the .pro file) is empty, then 
		for (count_t x=0; x<pro.size(); ++x) ind.push_back(x);  //put every clone in the vector ind.
	};

	std::vector <float_t> sum_gofs(ind.size() );
	std::vector <float_t> gofs_read(ind.size() );
	models model;
	allele_stat buffer_mle[BUFFER_SIZE]; 
	Locus buffer_site[BUFFER_SIZE];
	uint32_t all_read=0;
	size_t line_id = 0;
	while (true) {			//reads the next line of the pro file. pro.read() retuerns 0
		uint32_t c = 0, readed = 0;
		if (mpi_select(line_id, taskid, numtasks))
		{
			for (uint32_t x = 0; x < BUFFER_SIZE; ++x) {
				c = readed;				//Turn on the ability to read data from all clones in 
				if (pro.read(buffer_site[c]) != EOF) {
					readed++;	//reads the next line of the pro file. pro.read() retuerns 0
					line_id++;
					do_estimate(&buffer_mle[c], buffer_site[c], model, ind, sum_gofs, gofs_read, MIN, EMLMIN, MINGOF, MAXPITCH);
				}
			}
			if (taskid > 0)
			{
				std::cerr << "Sending from task " << taskid << " - last line " << line_id << '\n';
				std::ostringstream ost;
				write(ost, readed, pro, pro_out, buffer_mle, buffer_site, MINGOF, A);
				MPI_Send((void *)ost.str().c_str(), ost.str().length(), MPI_CHAR, 0, MAPGD_LINES, MPI_COMM_WORLD);
			}
			else
			{
				std::cerr << "Task 0 writing - last line " << line_id << '\n';
				write(*out, readed, pro, pro_out, buffer_mle, buffer_site, MINGOF, A);
				for (int i = 1; i < numtasks; ++i)
				{
					std::cerr << "Recieving from task " << i << "\n";
					std::istringstream ist(mpi_recieve_string(i));
					*out << ist.str();
				}
			}
		}
		else
		{
			// skip lines until our next stride
			for (uint32_t x = 0; x < BUFFER_SIZE; ++x) {
				c = readed;	
				if (pro.read(buffer_site[c]) != EOF) {
					readed++;
					line_id++;
				}
			}
		}
		if (readed!=BUFFER_SIZE)
		{
			if (taskid > 0)
			{
				std::ostringstream ost;
				write(ost, readed, pro, pro_out, buffer_mle, buffer_site, MINGOF, A);
				MPI_Send("", 0, MPI_CHAR, 0, MAPGD_LINES, MPI_COMM_WORLD);
			}
			break;
		}
		all_read+=readed;
		if (all_read>stop){break;}
	}
	if (taskid == 0)
	{
		// accumulate sum_gofs and gofs_read from all the running tasks
		for (int i = 1; i < numtasks; ++i)
		{
			std::vector<float_t> recv(ind.size() * 2);
			MPI_Recv(recv.data(), recv.size(), MPI_FLOAT, MPI_ANY_SOURCE, MAPGD_GOFS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			std::cerr << "Recieved GOF data\n";
			for (size_t x = 0; x < ind.size(); ++x)
			{
				gofs_read[x] += recv[x];
				sum_gofs[x] += recv[x + ind.size()];
			}
		}
		for (size_t x = 0; x<ind.size(); ++x)
			*out << "@" << pro.get_sample_name(ind[x]) << ":" << sum_gofs[x] / (float_t(gofs_read[x])) << std::endl;
	}
	else
	{
		std::vector<float_t> send(ind.size()*2);
		std::copy(gofs_read.begin(), gofs_read.end(), send.begin());
		std::copy(sum_gofs.begin(), sum_gofs.end(), send.begin()+ind.size());
		std::cerr << "Task " << taskid << " sending GOF data\n";
		MPI_Send(send.data(), send.size(), MPI_FLOAT, 0, MAPGD_GOFS, MPI_COMM_WORLD);

	}
	pro.close();
	if (outFile.is_open()) outFile.close();		//Closes outFile iff outFile is open.
	if (pro_out.is_open()) pro_out.close();		//Closes pro_out iff pro_out is open.

	MPI_Finalize();

	return 0;					//Since everything worked, return 0!.
}
