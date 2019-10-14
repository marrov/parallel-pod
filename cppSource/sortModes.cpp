#include <cstring>
#include <ctime>
#include <Eigen/Dense>
#include "ezOptionParser.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <string>
#include <sys/stat.h>
#include <sys/time.h>
#include <vector>

using namespace Eigen;

void sortrows(std::vector<std::vector<double>> &matrix, int col)
{
    std::sort(matrix.begin(),
              matrix.end(),
              [col](const std::vector<double> &lhs, const std::vector<double> &rhs) {
                  return lhs[col] > rhs[col];
              });
}

void pod(ez::ezOptionParser &opt)
{
    double start, end;

    std::cout << "Starting reconstruction routines \n"
              << std::endl;

    // Rows of matrix (number of points)
    long long MSIZE;
    opt.get("-p")->getLongLong(MSIZE);

    // Size of the variable (1 if scalar, 3 if vector)
    long long VSIZE;
    opt.get("-v")->getLongLong(VSIZE);

    // Size of input POD modes (number of modes to use in reconstruction)
    long long NSIZE;
    opt.get("-nm")->getLongLong(NSIZE);

    // Number of parallel threads
    long long PSIZE;
    opt.get("-np")->getLongLong(PSIZE);

    std::string dir_sort;
    opt.get("--sort")->getString(dir_sort);

    std::string dir_mode;
    opt.get("--mode")->getString(dir_mode);

    std::string dir_ref;
    opt.get("--reference")->getString(dir_ref);

    // READING MODE FILES

    omp_set_num_threads(PSIZE);
    start = omp_get_wtime();
    std::cout << "Reading modes..." << std::flush;
    MatrixXd m = MatrixXd::Zero(MSIZE, NSIZE * VSIZE);
    std::string xyz = "xyz_";
#pragma omp parallel
#pragma omp for
    for (size_t k = 0; k < NSIZE; k++)
    {
        for (size_t j = 0; j < VSIZE; j++)
        {
            std::ifstream readMode(dir_mode + "/mode_U" + xyz.at(j) + xyz.at(3) + std::to_string(k) + ".dat");

            if (readMode.is_open())
            {
                for (size_t i = 0; i < MSIZE; i++)
                {
                    readMode >> m(i, j + (VSIZE * k));
                }
            }
            else
            {
                std::cout << "Unable to open file" << std::endl;
            }
            readMode.close();
        }
    }
    end = omp_get_wtime();
    std::cout << "\t\t\t\t Done in " << end - start << "s \n"
              << std::endl;

    // READING REFERENCE FILE

    start = omp_get_wtime();
    std::cout << "Reading reference coordinate points..." << std::endl;
    MatrixXd r = MatrixXd::Zero(MSIZE, VSIZE);
    std::ifstream readRef(dir_ref + "pointCloud_reference.xy");
    double dummy = 0;

    if (readRef.is_open())
    {
        for (size_t i = 0; i < MSIZE; i++)
        {
            for (size_t j = 0; j < VSIZE + 1; j++)
            {
                if (j >= VSIZE)
                {
                    readRef >> dummy;
                }
                else
                {
                    readRef >> r(i, j);
                }
            }
        }
    }
    else
    {
        std::cout << "Unable to open file" << std::endl;
    }

    readRef.close();

    std::cout << r << std::endl;

    //sortrows(r,2);

    // SORT ALGORITHM
}

void Usage(ez::ezOptionParser &opt)
{
    std::string usage;
    opt.getUsage(usage);
    std::cout << usage;
};

int main(int argc, const char *argv[])
{
    ez::ezOptionParser opt;

    opt.overview = "POD routine";
    opt.syntax = "Perform POD (Proper Orthogonal Decomposiztion) using [INPUTS] ...";
    opt.example = "Add example \n\n";
    opt.footer = "POD routine. \nThis program is free and without warranty.\n";

    opt.add(
        "",                            // Default.
        0,                             // Required?
        0,                             // Number of args expected.
        0,                             // Delimiter if expecting multiple args.
        "Display usage instructions.", // Help description.
        "-h",                          // Flag token.
        "-help",                       // Flag token.
        "--help",                      // Flag token.
        "--usage"                      // Flag token.
    );

    opt.add(
        "",                  // Default.
        1,                   // Required?
        1,                   // Number of args expected.
        0,                   // Delimiter if expecting multiple args.
        "Output directory.", // Help description.
        "-s",                // Flag token.
        "-srt",              // Flag token.
        "-sort",             // Flag token.
        "--sort"             // Flag token.
    );

    opt.add(
        "",                 // Default.
        1,                  // Required?
        1,                  // Number of args expected.
        0,                  // Delimiter if expecting multiple args.
        "Modes directory.", // Help description.
        "-m",               // Flag token.
        "-md",              // Flag token.
        "-mode",            // Flag token.
        "--mode"            // Flag token.
    );

    opt.add(
        "",                                       // Default.
        1,                                        // Required?
        1,                                        // Number of args expected.
        0,                                        // Delimiter if expecting multiple args.
        "Point coorfinates reference directory.", // Help description.
        "-r",                                     // Flag token.
        "-ref",                                   // Flag token.
        "-reference",                             // Flag token.
        "--reference"                             // Flag token.
    );

    // Validator for unisgned long long (8 bytes).
    ez::ezOptionValidator *vU8 = new ez::ezOptionValidator("u8");
    opt.add(
        "",                               // Default.
        1,                                // Required?
        1,                                // Number of args expected.
        0,                                // Delimiter if expecting multiple args.
        "Number of points per snapshot.", // Help description.
        "-p",                             // Flag token.
        vU8                               //Validate input
    );

    opt.add(
        "",                            // Default.
        1,                             // Required?
        1,                             // Number of args expected.
        0,                             // Delimiter if expecting multiple args.
        "Number of values per point.", // Help description.
        "-v",                          // Flag token.
        vU8                            //Validate input
    );

    opt.add(
        "",                        // Default.
        1,                         // Required?
        1,                         // Number of args expected.
        0,                         // Delimiter if expecting multiple args.
        "Number of modes to read", // Help description.
        "-nm",                     // Flag token.
        vU8                        //Validate input
    );

    opt.add(
        "",                            // Default.
        1,                             // Required?
        1,                             // Number of args expected.
        0,                             // Delimiter if expecting multiple args.
        "Number of parallel threads.", // Help description.
        "-np",                         // Flag token.
        vU8                            //Validate input
    );

    // Perform the actual parsing of the command line.
    opt.parse(argc, argv);

    if (opt.isSet("-h"))
    {
        Usage(opt);
        return 1;
    }

    // Perform validations of input parameters.
    //
    // Check if directories exist.
    std::array<std::string, 3> dirflags = {"-s", "-m"};
    for (auto &dirflag : dirflags)
    {
        if (opt.isSet(dirflag))
        {
            std::string inputdir;
            struct stat info;
            opt.get(dirflag.c_str())->getString(inputdir);

            if (stat(inputdir.c_str(), &info) != 0)
            {
                std::cerr << "ERROR: " << inputdir << " does not exist.\n\n";
                return 1;
            }

            if (!(info.st_mode & S_IFDIR))
            { // S_ISDIR() doesn't exist on my windows
                std::cerr << "ERROR: " << inputdir << " is not a directory.\n\n";
                return 1;
            }
        }
    }

    std::vector<std::string> badOptions;
    int i;
    if (!opt.gotRequired(badOptions))
    {
        for (i = 0; i < badOptions.size(); ++i)
            std::cerr << "ERROR: Missing required option " << badOptions[i] << ".\n\n";

        Usage(opt);
        return 1;
    }

    if (!opt.gotExpected(badOptions))
    {
        for (i = 0; i < badOptions.size(); ++i)
            std::cerr << "ERROR: Got unexpected number of arguments for option " << badOptions[i] << ".\n\n";

        Usage(opt);
        return 1;
    }

    std::string firstArg;
    if (opt.firstArgs.size() > 0)
        firstArg = *opt.firstArgs[0];

    pod(opt);
}