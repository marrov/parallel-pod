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
#include <iterator>
#include <algorithm>

using namespace Eigen;

/*
Struct to hold information about the point cloud files.
*/
struct pointCloudFileInfo
{
    long rows;
    long columns;
};

/*
Function for splitting strings into a vector datatype.
*/
std::vector<std::string> split_string(const std::string &s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

/*
Function for reading the input file containing all the itme entries to use.
*/
std::vector<std::string> read_timefile(const std::string tfile)
{
    /* Open and read file with list of times */
    std::ifstream timefile(tfile);
    std::vector<std::string> t_entr;

    if (timefile.is_open())
    {
        /* Count number of words in file. This should correspond to the
        number of snapshots to use. */
        std::string line;

        while (getline(timefile, line))
        {
            t_entr.push_back(line);
        }
    }
    timefile.close();

    return t_entr;
}

/*
Parse the point cloud files and populate matrix with data.
*/
pointCloudFileInfo read_pcfs_to_matrix(MatrixXd *m,
                                       const std::vector<std::string> *fvec,
                                       const long no_cols,
                                       //std::vector<int> cols_to_read)
                                       const long offset)
{
    pointCloudFileInfo pointCloudRefFileInfo;
    bool verbose = false;

    if (verbose)
    {
        std::ostream_iterator<std::string> out_it(std::cout, "\n");
        std::copy(fvec->begin(), fvec->end(), out_it);
    }

    /* Establish a reference number of points for checking the problem size.
    The size is determined from the point cloud file in the first time directory. */
    pointCloudRefFileInfo.rows = 0;
    std::string ref_fname = *(fvec->begin());

    std::ifstream ref_file(ref_fname);
    if (ref_file.is_open())
    {
        std::string unused_line;

        /* Get a line and count the number of columns by splitting at spaces ' '.
        Increment the row counter. */
        getline(ref_file, unused_line);
        pointCloudRefFileInfo.rows++;
        auto tokens = split_string(unused_line, ' ');
        pointCloudRefFileInfo.columns = tokens.size();

        /* Count the rest of the rows. */
        while (getline(ref_file, unused_line))
            pointCloudRefFileInfo.rows++;
    }
    ref_file.close();

    if (verbose)
    {
        std::cout << "From file " << ref_fname << std::endl;
        std::cout << "found " << pointCloudRefFileInfo.rows << " points" << std::endl;
    }

    /* Number of time samples to consider is determined from the length of the
    vector of files */
    long TSIZE = fvec->size();

    /* Define matrix to store the file content */
    *m = MatrixXd::Zero(pointCloudRefFileInfo.rows * no_cols, TSIZE);
    double dummy = 0;
#pragma omp parallel
#pragma omp for
    for (size_t snapshot = 0; snapshot < TSIZE; snapshot++)
    {
        std::ifstream file((*fvec)[snapshot]);

        if (file.is_open())
        {
            for (size_t row_idx = 0; row_idx < pointCloudRefFileInfo.rows; row_idx++)
            {
                for (size_t col_no = 0, j = 0; col_no < (no_cols + offset); col_no++)
                {
                    if (col_no < offset)
                        file >> dummy;
                    else
                    {
                        file >> (*m)(row_idx + pointCloudRefFileInfo.rows * j, snapshot);
                        j++;
                    }
                }
            }

            /*
            std::string line;
            std::vector<std::string> tokens;
            char delimiter = '\n';
            unsigned long row_idx=0, col_idx=0;
            while (std::getline(file, line, delimiter))
            //for(;;std::getline(file, line, delimiter))
            {
                tokens = split_string(line, ' ');

                // for(std::vector<std::string>::iterator it = tokens.begin(); 
                //     it != tokens.end(); ++it)
                for(col_idx=0; col_idx<no_cols; col_idx++)
                {
                    auto tmp_val = tokens[(cols_to_read[col_idx]-1)];
                    (*m)(row_idx + pointCloudRefFileInfo.rows*col_idx, snapshot) = std::stod(tmp_val);
                }
            }
            */

            file.close();
        }
        else
        {
            std::cerr << "Unable to open file " << (*fvec)[snapshot] << std::endl;
        }
    }

    return pointCloudRefFileInfo;
}

void pod(ez::ezOptionParser &opt)
{
    double start, end;

    std::cout << "Starting POD routines \n"
              << std::endl;

    // Size of the variable (1 if scalar, 3 if vector)
    long long VSIZE;
    opt.get("-v")->getLongLong(VSIZE);

    // Offset to (column number of) the first column to read from.
    long long OFFSET;
    opt.get("-co")->getLongLong(OFFSET);

    // Size of input POD modes (number of modes to write)
    long long NSIZE;
    opt.get("-nm")->getLongLong(NSIZE);

    // Number of parallel threads
    long long PSIZE;
    opt.get("-np")->getLongLong(PSIZE);
    omp_set_num_threads(PSIZE);

    std::string dir_input;
    opt.get("--input")->getString(dir_input);

    std::string dir_chronos;
    opt.get("--chronos")->getString(dir_chronos);

    std::string dir_mode;
    opt.get("--mode")->getString(dir_mode);

    std::string tname;
    opt.get("-tf")->getString(tname);

    std::string pcfname;
    opt.get("-pcfn")->getString(pcfname);

    // GENERATING TIME STRING
    std::vector<std::string> t;

    t = read_timefile(tname);
    long TSIZE = t.size(); // Determine number of snapshots from time entry list

    /* Check whether the requested number of modes to write is larger than
    the number of snapshots. If so, set the number of modes to the number of
    snapshots. */
    if (NSIZE > TSIZE)
    {
        std::cout << "Modes to write exceed available snapshots. Adjusted.\n"
                  << std::flush;
        NSIZE = TSIZE;
    }

    // READING INPUT FILES

    /* Build a string array with all file names to be read into the matrix */
    std::vector<std::string> pcfs;
    for (std::vector<std::string>::iterator it = t.begin(); it != t.end(); ++it)
    {
        std::string name_temp = dir_input + "/" + *it + "/" + pcfname;
        pcfs.push_back(name_temp);
    }

    MatrixXd m;
    start = omp_get_wtime();
    std::cout << "Reading files..." << std::flush;
    // auto REF_MSIZE = read_pcfs_to_matrix(&m, &pcfs, (long)VSIZE);
    // std::vector<int> cols_to_read = {1,2,3};
    // auto REF_MSIZE = read_pcfs_to_matrix(&m, &pcfs, (long)VSIZE, cols_to_read);
    //int offset = 0;
    //auto REF_MSIZE = read_pcfs_to_matrix(&m, &pcfs, (long)VSIZE, (long)OFFSET);
    auto pointCloudInfo = read_pcfs_to_matrix(&m, &pcfs, (long)VSIZE, (long)OFFSET);
    auto REF_MSIZE = pointCloudInfo.rows;
    end = omp_get_wtime();
    std::cout << "\t\t\t\t Done in " << end - start << "s \n"
              << std::endl;

    std::cout << "File contains " << pointCloudInfo.rows << " rows and " << pointCloudInfo.columns << " columns. "
              << "Read data from columns " << (OFFSET + 1) << " to " << (OFFSET + VSIZE) << ".\n"
              << std::endl;

    // COMPUTING NORMALISED PROJECTION MATRIX

    start = omp_get_wtime();
    std::cout << "Computing projection matrix..." << std::flush;
    MatrixXd pm = MatrixXd::Zero(TSIZE, TSIZE);
    pm = (1.0 / TSIZE) * m.transpose() * m;
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
                    spm(i, j) += g(idx) * pmExt(i + k + SPOD_Nf, j + k + SPOD_Nf);
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
    MatrixXd podx = MatrixXd::Zero(REF_MSIZE, TSIZE);
    MatrixXd pody = MatrixXd::Zero(REF_MSIZE, TSIZE);
    MatrixXd podz = MatrixXd::Zero(REF_MSIZE, TSIZE);
#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < TSIZE; i++)
    {
        for (size_t j = 0; j < TSIZE; j++)
        {
            //podx.col(i) = podx.col(i) + (1.0 / (eigval(i) * TSIZE)) * sqrt(eigval(i) * TSIZE) * eigvec(j, i) * m.block(0 * MSIZE, j, MSIZE, 1);
            //pody.col(i) = pody.col(i) + (1.0 / (eigval(i) * TSIZE)) * sqrt(eigval(i) * TSIZE) * eigvec(j, i) * m.block(1 * MSIZE, j, MSIZE, 1);
            //podz.col(i) = podz.col(i) + (1.0 / (eigval(i) * TSIZE)) * sqrt(eigval(i) * TSIZE) * eigvec(j, i) * m.block(2 * MSIZE, j, MSIZE, 1);

            podx.col(i) += (eigvec(j, i) / sqrt(eigval(i) * TSIZE)) * m.block(0 * REF_MSIZE, j, REF_MSIZE, 1);
            pody.col(i) += (eigvec(j, i) / sqrt(eigval(i) * TSIZE)) * m.block(1 * REF_MSIZE, j, REF_MSIZE, 1);
            podz.col(i) += (eigvec(j, i) / sqrt(eigval(i) * TSIZE)) * m.block(2 * REF_MSIZE, j, REF_MSIZE, 1);
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
        "",                             // Default.
        1,                              // Required?
        1,                              // Number of args expected.
        0,                              // Delimiter if expecting multiple args.
        "File with time snapshot list", // Help description.
        "-tf"                           // Flag token.
    );

    opt.add(
        "",                                             // Default.
        1,                                              // Required?
        1,                                              // Number of args expected.
        0,                                              // Delimiter if expecting multiple args.
        "Point cloud file name (in time directories).", // Help description.
        "-pcfn"                                         // Flag token.
    );

    opt.add(
        "",                                                                // Default.
        1,                                                                 // Required?
        1,                                                                 // Number of args expected.
        0,                                                                 // Delimiter if expecting multiple args.
        "Directory where the time directories reside as sub-directories.", // Help description.
        "-i",                                                              // Flag token.
        "-inp",                                                            // Flag token.
        "-input",                                                          // Flag token.
        "--input"                                                          // Flag token.
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
        "",                            // Default.
        1,                             // Required?
        1,                             // Number of args expected.
        0,                             // Delimiter if expecting multiple args.
        "Number of values per point.", // Help description.
        "-v",                          // Flag token.
        vU8                            // Validate input
    );

    opt.add(
        "0",                                                 // Default.
        0,                                                   // Required?
        1,                                                   // Number of args expected.
        0,                                                   // Delimiter if expecting multiple args.
        "Point cloud file column offset (for reading data)", // Help description.
        "-co",                                               // Flag token.
        vU8                                                  // Validate input
    );

    opt.add(
        "",                          // Default.
        1,                           // Required?
        1,                           // Number of args expected.
        0,                           // Delimiter if expecting multiple args.
        "Number of modes to write.", // Help description.
        "-nm",                       // Flag token.
        vU8                          // Validate input
    );

    opt.add(
        "",                            // Default.
        1,                             // Required?
        1,                             // Number of args expected.
        0,                             // Delimiter if expecting multiple args.
        "Number of parallel threads.", // Help description.
        "-np",                         // Flag token.
        vU8                            // Validate input
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
    return 0;
}
