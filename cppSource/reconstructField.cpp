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

void pod(ez::ezOptionParser &opt)
{
    double start, end;

    std::cout << "Starting reconstruction routines \n"
              << std::endl;

    // Rows of matrix (number of points)
    long long MSIZE;
    opt.get("-p")->getLongLong(MSIZE);

    // Size of time data (number of snapshots)
    long long TSIZE;
    opt.get("-s")->getLongLong(TSIZE);

    // Size of the variable (1 if scalar, 3 if vector)
    long long VSIZE;
    opt.get("-v")->getLongLong(VSIZE);

    // Size of input POD modes (number of modes to use in reconstruction)
    long long NSIZE;
    opt.get("-nm")->getLongLong(NSIZE);

    // Number of parallel threads
    long long PSIZE;
    opt.get("-np")->getLongLong(PSIZE);

    // Number of reconstructed time instants to output
    long long RSIZE;
    opt.get("-nr")->getLongLong(RSIZE);

    std::string dir_rec;
    opt.get("--reconstruct")->getString(dir_rec);

    std::string dir_chronos;
    opt.get("--chronos")->getString(dir_chronos);

    std::string dir_mode;
    opt.get("--mode")->getString(dir_mode);

    // READING MODE FILES

    omp_set_num_threads(PSIZE);
    start = omp_get_wtime();
    std::cout << "Reading modes..." << std::flush;
    MatrixXd m = MatrixXd::Zero(MSIZE * VSIZE, NSIZE);
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
                    readMode >> m(i + MSIZE * j, k);
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

    // READING CHRONOS FILES
    start = omp_get_wtime();
    std::cout << "Reading chronos..." << std::flush;
    MatrixXd c = MatrixXd::Zero(TSIZE, NSIZE);
#pragma omp parallel
#pragma omp for
    for (size_t j = 0; j < NSIZE; j++)
    {

        std::ifstream readChronos(dir_chronos + "/chronos_" + std::to_string(j) + ".dat");

        if (readChronos.is_open())
        {
            for (size_t i = 0; i < TSIZE; i++)
            {
                readChronos >> c(i, j);
            }
        }
        else
        {
            std::cout << "Unable to open file" << std::endl;
        }
        readChronos.close();
    }
    end = omp_get_wtime();
    std::cout << "\t\t\t\t Done in " << end - start << "s \n"
              << std::endl;

    // COMPUTE RECONSTRUCTED FIELDS

    start = omp_get_wtime();
    std::cout << "Computing reconstructed fields..." << std::flush;
    MatrixXd recx = MatrixXd::Zero(MSIZE, RSIZE);
    MatrixXd recy = MatrixXd::Zero(MSIZE, RSIZE);
    MatrixXd recz = MatrixXd::Zero(MSIZE, RSIZE);
#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < RSIZE; i++)
    {
        for (size_t j = 0; j < NSIZE; j++)
        {
            recx.col(i) = recx.col(i) + c(i, j) * m.block(0 * MSIZE, j, MSIZE, 1);
            recy.col(i) = recy.col(i) + c(i, j) * m.block(1 * MSIZE, j, MSIZE, 1);
            recz.col(i) = recz.col(i) + c(i, j) * m.block(2 * MSIZE, j, MSIZE, 1);
        }
    }
    end = omp_get_wtime();
    std::cout << "\t\t Done in " << end - start << "s \n"
              << std::endl;

    // WRITE RECONSTRUCTED FIELDS

    start = omp_get_wtime();
    std::cout << "Writing reconstructed fields..." << std::flush;
#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < RSIZE; i++)
    {
        for (size_t j = 0; j < VSIZE; j++)
        {
            std::ofstream writeField(dir_rec + "/U" + xyz.at(j) + xyz.at(3) + std::to_string(i) + ".dat");
            if (writeField.is_open())
            {
                if (j == 0)
                {
                    writeField << std::scientific << std::setprecision(10) << recx.col(i);
                }
                else if (j == 1)
                {
                    writeField << std::scientific << std::setprecision(10) << recy.col(i);
                }
                else if (j == 2)
                {
                    writeField << std::scientific << std::setprecision(10) << recz.col(i);
                }
            }
            writeField.close();
        }
    }
    end = omp_get_wtime();
    std::cout << "\t\t\t Done in " << end - start << "s \n"
              << std::endl;
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
        "",                 // Default.
        1,                  // Required?
        1,                  // Number of args expected.
        0,                  // Delimiter if expecting multiple args.
        "Input directory.", // Help description.
        "-r",               // Flag token.
        "-rec",             // Flag token.
        "-reconstruct",     // Flag token.
        "--reconstruct"     // Flag token.
    );

    opt.add(
        "",                   // Default.
        1,                    // Required?
        1,                    // Number of args expected.
        0,                    // Delimiter if expecting multiple args.
        "Chronos directory.", // Help description.
        "-c",                 // Flag token.
        "-chr",               // Flag token.
        "-chronos",           // Flag token.
        "--chronos"           // Flag token.
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
        "",                                       // Default.
        1,                                        // Required?
        1,                                        // Number of args expected.
        0,                                        // Delimiter if expecting multiple args.
        "Number of sequential snapshots to use.", // Help description.
        "-s",                                     // Flag token.
        vU8                                       //Validate input
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
        "",                                                                        // Default.
        1,                                                                         // Required?
        1,                                                                         // Number of args expected.
        0,                                                                         // Delimiter if expecting multiple args.
        "Number of modes to read (i.e. number of modes to use in reconstruction)", // Help description.
        "-nm",                                                                     // Flag token.
        vU8                                                                        //Validate input
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

    opt.add(
        "",                                                 // Default.
        1,                                                  // Required?
        1,                                                  // Number of args expected.
        0,                                                  // Delimiter if expecting multiple args.
        "Number of reconstructed time instants to output.", // Help description.
        "-nr",                                              // Flag token.
        vU8                                                 //Validate input
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
    std::array<std::string, 3> dirflags = {"-r", "-c", "-m"};
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

    // Check if number of snapshots compared to modes to read.
    // Size of time data (number of snapshots)
    long long TSIZE;
    opt.get("-s")->getLongLong(TSIZE);
    // Size of output POD modes (number of modes to read)
    long long NSIZE;
    opt.get("-nm")->getLongLong(NSIZE);
    if (TSIZE <= NSIZE)
    {
        std::cerr << "ERROR: Number of modes to read must be less or equal to number of snapshots used.\n\n";

        Usage(opt);
        return 1;
    }

    long long RSIZE;
    opt.get("-nm")->getLongLong(RSIZE);
    if (TSIZE <= RSIZE)
    {
        std::cerr << "ERROR: Number of reconstructed time instants to read must be less or equal to number of snapshots used.\n\n";

        Usage(opt);
        return 1;
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