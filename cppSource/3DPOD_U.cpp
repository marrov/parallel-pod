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

    std::cout << "Starting POD routines \n"
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

    // Size of input POD modes (number of modes to write)
    long long NSIZE;
    opt.get("-nm")->getLongLong(NSIZE);

    // Number of parallel threads
    long long PSIZE;
    opt.get("-np")->getLongLong(PSIZE);

    std::string dir_input;
    opt.get("--input")->getString(dir_input);

    std::string dir_chronos;
    opt.get("--chronos")->getString(dir_chronos);

    std::string dir_mode;
    opt.get("--mode")->getString(dir_mode);

    // GENERATING TIME STRING

    std::vector<std::string> t(TSIZE);
    std::string tname = "times_pointCloud.txt";
    std::ifstream timefile(dir_input + "/" + tname);

    for (size_t i = 0; i < TSIZE; i++)
    {
        if (timefile.is_open())
        {
            timefile >> t[i];
            while (t[i].back() == '0') // Remove trailing zeros
            {
                t[i].pop_back();
            }
        }
    }

    // READING INPUT FILES

    omp_set_num_threads(PSIZE);
    start = omp_get_wtime();
    std::cout << "Reading files..." << std::flush;
    MatrixXd m = MatrixXd::Zero(MSIZE * VSIZE, TSIZE);
#pragma omp parallel
#pragma omp for
    for (size_t k = 0; k < TSIZE; k++)
    {
        std::string dir = dir_input + "/pointCloud/" + t[k] + "/pointCloud_U.xy";
        std::ifstream file(dir);

        if (file.is_open())
        {
            for (size_t i = 0; i < MSIZE; i++)
            {
                for (size_t j = 0; j < VSIZE; j++)
                {
                    file >> m(i + MSIZE * j, k);
                }
            }

            // IF VERBOSE do print which file
            //std::cout << "Finished reading file " + std::to_string(k + 1) + " \t of " + std::to_string(TSIZE) << " by thread " << omp_get_thread_num() << std::endl;

            file.close();
        }
        else
        {
            std::cout << "Unable to open file" << std::endl;
        }
    }
    end = omp_get_wtime();
    std::cout << "\t\t\t\t Done in " << end - start << "s \n"
              << std::endl;

    // COMPUTING NORMALISED PROJECTION MATRIX

    start = omp_get_wtime();
    std::cout << "Computing projection matrix..." << std::flush;
    MatrixXd pm = MatrixXd::Zero(TSIZE, TSIZE);
#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < TSIZE; i++)
    {
        for (size_t j = 0; j < TSIZE; j++)
        {
            pm(i, j) = (1.0 / TSIZE) * (m.col(i).dot(m.col(j).transpose()));
        }
    }
    end = omp_get_wtime();
    std::cout << "\t\t\t Done in " << end - start << "s \n"
              << std::endl;

    // APPLY SPECTRAL POD FILTER IF DESIRED

    int SPOD_Fl = 0; // Flag to apply SPOD - 1 = on
    int SPOD_Ty = 2; // Filter type - 1 = box, 2 = gauss
    int SPOD_Nf = 5; // Filter size, Nf

    if (SPOD_Fl == 1)
    {
        start = omp_get_wtime();
        std::cout << "Filtering projection matrix for SPOD..." << std::flush;

        int nfSize = 2 * SPOD_Nf + 1;

        VectorXd g = VectorXd::Ones(nfSize);

        if (SPOD_Ty == 2)
        {
            VectorXd gauss = VectorXd::LinSpaced(nfSize, -2.285, 2.285);
            g = exp(-square(gauss.array()));
        }

        g = g / g.sum();

        size_t idx = 0;

        MatrixXd spm = MatrixXd::Zero(TSIZE, TSIZE);
        MatrixXd pmExt = MatrixXd::Zero(TSIZE * 3, TSIZE * 3);

        pmExt = pm.replicate(3, 3).block(TSIZE - SPOD_Nf, TSIZE - SPOD_Nf, TSIZE + 2 * SPOD_Nf, TSIZE + 2 * SPOD_Nf);

        for (int i = 0; i < TSIZE; i++)
        {
            for (int j = 0; j < TSIZE; j++)
            {
                for (int k = -SPOD_Nf; k < SPOD_Nf + 1; k++)
                {
                    spm(i, j) = spm(i, j) + g(idx) * pmExt(i + k + SPOD_Nf, j + k + SPOD_Nf);
                    idx++;
                }
                idx = 0;
            }
        }

        pm = spm;

        end = omp_get_wtime();
        std::cout << "\t\t Done in " << end - start << "s \n"
                  << std::endl;
    }

    // COMPUTING SORTED EIGENVALUES AND EIGENVECTORS

    start = omp_get_wtime();
    std::cout << "Computing eigenvalues and eigenvectors..." << std::flush;
    SelfAdjointEigenSolver<MatrixXd> eigensolver(pm);
    if (eigensolver.info() != Success)
        abort();

    VectorXd eigval = eigensolver.eigenvalues().reverse();
    MatrixXd eigvec = eigensolver.eigenvectors().rowwise().reverse();
    end = omp_get_wtime();
    std::cout << "\t Done in " << end - start << "s \n"
              << std::endl;

    // COMPUTING POD MODES

    start = omp_get_wtime();
    std::cout << "Computing POD modes..." << std::flush;
    MatrixXd podx = MatrixXd::Zero(MSIZE, TSIZE);
    MatrixXd pody = MatrixXd::Zero(MSIZE, TSIZE);
    MatrixXd podz = MatrixXd::Zero(MSIZE, TSIZE);
#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < TSIZE; i++)
    {
        for (size_t j = 0; j < TSIZE; j++)
        {
            podx.col(i) = podx.col(i) + (1.0 / (eigval(i) * TSIZE)) * sqrt(eigval(i) * TSIZE) * eigvec(j, i) * m.block(0 * MSIZE, j, MSIZE, 1);
            pody.col(i) = pody.col(i) + (1.0 / (eigval(i) * TSIZE)) * sqrt(eigval(i) * TSIZE) * eigvec(j, i) * m.block(1 * MSIZE, j, MSIZE, 1);
            podz.col(i) = podz.col(i) + (1.0 / (eigval(i) * TSIZE)) * sqrt(eigval(i) * TSIZE) * eigvec(j, i) * m.block(2 * MSIZE, j, MSIZE, 1);
        }
    }
    end = omp_get_wtime();
    std::cout << "\t\t\t\t Done in " << end - start << "s \n"
              << std::endl;

    // WRITING SORTED EIGENVALUES

    start = omp_get_wtime();
    std::cout << "Writing eigenvalues..." << std::flush;
    std::ofstream writeEigval(dir_chronos + "/A.txt");
    if (writeEigval.is_open())
    {
        writeEigval << std::scientific << std::setprecision(10) << eigval;
        writeEigval.close();
    }
    end = omp_get_wtime();
    std::cout << "\t\t\t\t Done in " << end - start << "s \n"
              << std::endl;

    // WRITING CHRONOS

    start = omp_get_wtime();
    std::cout << "Writing chronos..." << std::flush;
#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < NSIZE; i++)
    {
        std::ofstream writeChronos(dir_chronos + "/chronos_" + std::to_string(i) + ".dat");
        for (size_t j = 0; j < TSIZE; j++)
        {
            writeChronos << std::scientific << std::setprecision(10) << sqrt(eigval(i) * TSIZE) * eigvec(j, i) << '\n';
        }
        writeChronos.close();
    }
    end = omp_get_wtime();
    std::cout << "\t\t\t\t Done in " << end - start << "s \n"
              << std::endl;

    // WRITING POD MODES

    start = omp_get_wtime();
    std::cout << "Writing POD modes..." << std::flush;
    std::string xyz = "xyz_";
#pragma omp parallel
#pragma omp for
    for (size_t j = 0; j < VSIZE; j++)
    {
        for (size_t i = 0; i < NSIZE; i++)
        {
            std::ofstream writeMode(dir_mode + "/mode_U" + xyz.at(j) + xyz.at(3) + std::to_string(i) + ".dat");
            if (writeMode.is_open())
            {
                if (j == 0)
                {
                    writeMode << std::scientific << std::setprecision(6) << podx.col(i);
                }
                else if (j == 1)
                {
                    writeMode << std::scientific << std::setprecision(6) << pody.col(i);
                }
                else if (j == 2)
                {
                    writeMode << std::scientific << std::setprecision(6) << podz.col(i);
                }
                writeMode.close();
            }
        }
    }
    end = omp_get_wtime();
    std::cout << "\t\t\t\t Done in " << end - start << "s \n"
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
        "-i",               // Flag token.
        "-inp",             // Flag token.
        "-input",           // Flag token.
        "--input"           // Flag token.
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
        "",                                                                             // Default.
        1,                                                                              // Required?
        1,                                                                              // Number of args expected.
        0,                                                                              // Delimiter if expecting multiple args.
        "Number of modes to write. Must be less or equal to number of snapshots used.", // Help description.
        "-nm",                                                                          // Flag token.
        vU8                                                                             //Validate input
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
    std::array<std::string, 3> dirflags = {"-i", "-c", "-m"};
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

    // Check if number of snapshots compared to modes to write.
    // Size of time data (number of snapshots)
    long long TSIZE;
    opt.get("-s")->getLongLong(TSIZE);
    // Size of output POD modes (number of modes to write)
    long long NSIZE;
    opt.get("-nm")->getLongLong(NSIZE);
    if (TSIZE < NSIZE)
    {
        std::cerr << "ERROR: Number of modes to write must be less or equal to number of snapshots used.\n\n";

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