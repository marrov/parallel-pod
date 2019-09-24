#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <ctime>
#include "ezOptionParser.hpp"
#include <sys/stat.h>
#include <sys/stat.h>
//#include <array>

#define MSIZE 381600 // Rows of matrix (number of points)
#define TSIZE 5      // Size of time data (number of snapshots)
#define VSIZE 3      // Size of the variable (1 if scalar, 3 if vector)
#define NSIZE 5      // Size of output POD modes (number of modes to write)

using namespace Eigen;

void pod(ez::ezOptionParser& opt)
{
    std::cout << "Starting POD routines" << std::endl;

    clock_t start,end;

    MatrixXd m = MatrixXd::Zero(MSIZE * VSIZE, TSIZE);
    MatrixXd pm = MatrixXd::Zero(TSIZE, TSIZE);

    std::string dir_input;
    opt.get("--input")->getString(dir_input);

    std::string dir_chronos;
    opt.get("--chronos")->getString(dir_chronos);

    std::string dir_mode;
    opt.get("--mode")->getString(dir_mode);

    // READING INPUT FILES

    start = clock();
    for (size_t k = 0; k < TSIZE; k++)
    {
        //std::string dir = "input/U/";
        std::string fname = "U_" + std::to_string(k + 1) + ".dat";
        //std::ifstream file(dir + fname);
        std::ifstream file(dir_input + "/" + fname);

        if (file.is_open())
        {
            for (size_t i = 0; i < MSIZE; i++)
            {
                for (size_t j = 0; j < VSIZE; j++)
                {
                    file >> m(i + MSIZE * j, k);
                }
            }
            std::cout << "Finished reading " << fname << std::endl;

            file.close();
        }
        else
        {
            std::cout << " Unable to open file" << std::endl;
        }
    }
    end = clock();
    double runTime = (double)(end-start)/CLOCKS_PER_SEC;
    std::cout << "Reading files takes "<< runTime << "s"<< std::endl;


    // COMPUTING NORMALISED PROJECTION MATRIX

    start = clock();
    for (size_t i = 0; i < TSIZE; i++)
    {
        for (size_t j = 0; j < TSIZE; j++)
        {
            pm(i, j) = (1.0 / TSIZE) * (m.col(i).dot(m.col(j).transpose()));
        }
    }
    end = clock();
    runTime = (double)(end-start)/CLOCKS_PER_SEC;
    std::cout << "computing projection matrix takes "<< runTime << "s"<< std::endl;

    /*
    std::cout << "Normalised rojection Matrix:" << std::endl
              << std::setprecision(12) << pm << std::endl;
    */

    // COMPUTING SORTED EIGENVALUES AND EIGENVECTORS

    start = clock();
    SelfAdjointEigenSolver<MatrixXd> eigensolver(pm);
    if (eigensolver.info() != Success)
        abort();

    VectorXd eigval = eigensolver.eigenvalues().reverse();
    MatrixXd eigvec = eigensolver.eigenvectors().rowwise().reverse();
    end = clock();
    runTime = (double)(end-start)/CLOCKS_PER_SEC;
    std::cout << "eigen-decomposition takes "<< runTime << "s"<< std::endl;

    /*
    std::cout << "Sorted eigenvalues of Projection Matrix:\n"
              << eigval << std::endl;
    std::cout << "Matrix whose columns are the sorted eigenvectors of Projection Matrix:\n"
              << eigvec << std::endl;
    */

    // WRITING SORTED EIGENVALUES

    //std::ofstream writeEigval("chronos/A.txt");
    std::ofstream writeEigval(dir_chronos + "/A.txt");
    if (writeEigval.is_open())
    {
        writeEigval << std::scientific << std::setprecision(10) << eigval;
        writeEigval.close();
    }

    // COMPUTING POD MODES AND WRITING CHRONOS FOR EACH MODE

    MatrixXd podx = MatrixXd::Zero(MSIZE, TSIZE);
    MatrixXd pody = MatrixXd::Zero(MSIZE, TSIZE);
    MatrixXd podz = MatrixXd::Zero(MSIZE, TSIZE);

    start = clock();
    for (size_t i = 0; i < TSIZE; i++)
    {
        //std::string chronos = "chronos/chronos_" + std::to_string(i) + ".dat";
        std::string chronos = dir_chronos + "/chronos_" + std::to_string(i) + ".dat";
        std::ofstream writeChronos(chronos);

        for (size_t j = 0; j < TSIZE; j++)
        {
            if (i <= NSIZE)
            {
                writeChronos << std::scientific << std::setprecision(10) << sqrt(eigval(i) * TSIZE) * eigvec(j, i) << '\n';
            }

            podx.col(i) = podx.col(i) + (1.0 / (eigval(i) * TSIZE)) * sqrt(eigval(i) * TSIZE) * eigvec(j, i) * m.block<MSIZE, 1>(0 * MSIZE, j);
            pody.col(i) = pody.col(i) + (1.0 / (eigval(i) * TSIZE)) * sqrt(eigval(i) * TSIZE) * eigvec(j, i) * m.block<MSIZE, 1>(1 * MSIZE, j);
            podz.col(i) = podz.col(i) + (1.0 / (eigval(i) * TSIZE)) * sqrt(eigval(i) * TSIZE) * eigvec(j, i) * m.block<MSIZE, 1>(2 * MSIZE, j);
        }

        writeChronos.close();
    }
    end = clock();
    runTime = (double)(end-start)/CLOCKS_PER_SEC;
    std::cout << "computing POD modes takes "<< runTime << "s"<< std::endl;

    // WRITING POD MODES

    start = clock();
    for (size_t i = 0; i < NSIZE; i++)
    {
        //std::string modex = "mode/mode_Ux_" + std::to_string(i) + ".dat";
        std::string modex = dir_mode + "/mode_Ux_" + std::to_string(i) + ".dat";
        //std::string modey = "mode/mode_Uy_" + std::to_string(i) + ".dat";
        std::string modey = dir_mode + "/mode_Uy_" + std::to_string(i) + ".dat";
        //std::string modez = "mode/mode_Uz_" + std::to_string(i) + ".dat";
        std::string modez = dir_mode + "/mode_Uz_" + std::to_string(i) + ".dat";

        std::ofstream writeModex(modex);
        if (writeModex.is_open())
        {
            writeModex << std::scientific << std::setprecision(6) << podx.col(i);
            writeModex.close();
        }

        std::ofstream writeModey(modey);
        if (writeModey.is_open())
        {
            writeModey << std::scientific << std::setprecision(6) << pody.col(i);
            writeModey.close();
        }

        std::ofstream writeModez(modez);
        if (writeModez.is_open())
        {
            writeModez << std::scientific << std::setprecision(6) << podz.col(i);
            writeModez.close();
        }
    }
    end = clock();
    runTime = (double)(end-start)/CLOCKS_PER_SEC;
    std::cout << "writing POD modes takes "<< runTime << "s"<< std::endl;

}

void Usage(ez::ezOptionParser& opt) {
	std::string usage;
	opt.getUsage(usage);
	std::cout << usage;
};

int main(int argc, const char * argv[]){
    ez::ezOptionParser opt;

    opt.overview = "POD routine";
	opt.syntax = "Perform POD (Proper Orthogonal Decomposiztion) using [INPUTS] ...";
	opt.example = "Add example \n\n";
	opt.footer = "POD routine. \nThis program is free and without warranty.\n";

    opt.add(
		"", // Default.
		0, // Required?
		0, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Display usage instructions.", // Help description.
		"-h",     // Flag token. 
		"-help",  // Flag token.
		"--help", // Flag token.
		"--usage" // Flag token.
	);

    opt.add(
		"", // Default.
		1, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Input directory.", // Help description.
		"-i", // Flag token.
		"-inp", // Flag token.
		"-input", // Flag token.
		"--input" // Flag token.
	);

    opt.add(
		"", // Default.
		1, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Chronos directory.", // Help description.
		"-c", // Flag token.
		"-chr", // Flag token.
		"-chronos", // Flag token.
		"--chronos" // Flag token.
	);

    opt.add(
		"", // Default.
		1, // Required?
		1, // Number of args expected.
		0, // Delimiter if expecting multiple args.
		"Modes directory.", // Help description.
		"-m", // Flag token.
		"-md", // Flag token.
		"-mode", // Flag token.
		"--mode" // Flag token.
	);

    // Perform the actual parsing of the command line.
    opt.parse(argc, argv);

    if (opt.isSet("-h")) {
		Usage(opt);
		return 1;
	}

    // Perform validations of input parameters.
    //
    // Check if directories exist.
    std::array<std::string, 3> dirflags = {"-i", "-c", "-m"};
    for(auto& dirflag : dirflags) { 
        if (opt.isSet(dirflag)) {
	    	std::string inputdir;
            struct stat info;
	    	opt.get(dirflag.c_str())->getString(inputdir);

            if(stat(inputdir.c_str(), &info) !=0 ){
                std::cerr << "ERROR: " << inputdir << " does not exist.\n\n";
                return 1;
            }

            if(!(info.st_mode & S_IFDIR)){  // S_ISDIR() doesn't exist on my windows 
                std::cerr << "ERROR: " << inputdir << " is not a directory.\n\n";
                return 1;
            }
	    }
    }

    std::vector<std::string> badOptions;
	int i;
	if(!opt.gotRequired(badOptions)) {
		for(i=0; i < badOptions.size(); ++i)
			std::cerr << "ERROR: Missing required option " << badOptions[i] << ".\n\n";

		Usage(opt);
		return 1;
	}

	if(!opt.gotExpected(badOptions)) {
		for(i=0; i < badOptions.size(); ++i)
			std::cerr << "ERROR: Got unexpected number of arguments for option " << badOptions[i] << ".\n\n";
			
		Usage(opt);
		return 1;
	}

	std::string firstArg;
	if (opt.firstArgs.size() > 0)
		firstArg = *opt.firstArgs[0];

    pod(opt);
}