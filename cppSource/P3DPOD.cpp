#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <string>
#include <sys/time.h>
#include <vector>

#define MSIZE 381600 // Rows of matrix (number of points)
#define TSIZE 200    // Size of time data (number of snapshots)
#define VSIZE 3      // Size of the variable (1 if scalar, 3 if vector)
#define NSIZE 20     // Size of output POD modes (number of modes to write)
#define PSIZE 4      // Number of parallel threads

using namespace Eigen;

int main()
{
    // GENERATING TIME STRING

    std::vector<std::string> t(TSIZE);

    std::string timedir = "../testCase/input/times_pointCloud.txt";
    std::ifstream timefile(timedir);

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

    MatrixXd m = MatrixXd::Zero(MSIZE * VSIZE, TSIZE);

    omp_set_num_threads(PSIZE);

#pragma omp parallel
#pragma omp for
    for (size_t k = 0; k < TSIZE; k++)
    {

        std::string dir = "../testCase/input/pointCloud/" + t[k] + "/pointCloud_U.xy";
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
            std::cout << "Finished reading file " + std::to_string(k + 1) + " of \t" + std::to_string(TSIZE) << " by thread " << omp_get_thread_num() << std::endl;

            file.close();
        }
        else
        {
            std::cout << "Unable to open file" << std::endl;
        }
    }

    // COMPUTING NORMALISED PROJECTION MATRIX

    MatrixXd pm = MatrixXd::Zero(TSIZE, TSIZE);

    std::cout << "Computing projection matrix..." << std::flush;

#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < TSIZE; i++)
    {
        for (size_t j = 0; j < TSIZE; j++)
        {
            pm(i, j) = (1.0 / TSIZE) * (m.col(i).dot(m.col(j).transpose()));
        }
    }
    std::cout << " Done" << std::endl;

    // COMPUTING SORTED EIGENVALUES AND EIGENVECTORS

    std::cout << "Computing eigenvalues..." << std::flush;
    SelfAdjointEigenSolver<MatrixXd>
        eigensolver(pm);
    if (eigensolver.info() != Success)
        abort();

    VectorXd eigval = eigensolver.eigenvalues().reverse();
    MatrixXd eigvec = eigensolver.eigenvectors().rowwise().reverse();
    std::cout << " Done" << std::endl;

    // COMPUTING POD MODES

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
            podx.col(i) = podx.col(i) + (1.0 / (eigval(i) * TSIZE)) * sqrt(eigval(i) * TSIZE) * eigvec(j, i) * m.block<MSIZE, 1>(0 * MSIZE, j);
            pody.col(i) = pody.col(i) + (1.0 / (eigval(i) * TSIZE)) * sqrt(eigval(i) * TSIZE) * eigvec(j, i) * m.block<MSIZE, 1>(1 * MSIZE, j);
            podz.col(i) = podz.col(i) + (1.0 / (eigval(i) * TSIZE)) * sqrt(eigval(i) * TSIZE) * eigvec(j, i) * m.block<MSIZE, 1>(2 * MSIZE, j);
        }
    }
    std::cout << " Done" << std::endl;

    // WRITING SORTED EIGENVALUES

    std::cout << "Writing eigenvalues..." << std::flush;
    std::ofstream writeEigval("../testCase/output/chronos/A.dat");
    if (writeEigval.is_open())
    {
        writeEigval << std::scientific << std::setprecision(10) << eigval;
        writeEigval.close();
    }
    std::cout << " Done" << std::endl;

    // WRITING CHRONOS

    std::cout << "Writing chronos..." << std::flush;

#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < NSIZE; i++)
    {
        std::string chronos = "../testCase/output/chronos/chronos_" + std::to_string(i) + ".dat";
        std::ofstream writeChronos(chronos);
        for (size_t j = 0; j < TSIZE; j++)
        {
            writeChronos << std::scientific << std::setprecision(10) << sqrt(eigval(i) * TSIZE) * eigvec(j, i) << '\n';
        }
        writeChronos.close();
    }
    std::cout << " Done" << std::endl;

    // WRITING POD MODES

    std::cout << "Writing POD modes..." << std::flush;
    std::string xyz = "xyz_";
    std::string modeHead = "../testCase/output/mode/mode_U";

#pragma omp parallel
#pragma omp for
    for (size_t j = 0; j < VSIZE; j++)
    {
        for (size_t i = 0; i < NSIZE; i++)
        {
            std::string modeTail = std::to_string(i) + ".dat";

            std::ofstream writeMode(modeHead + xyz.at(j) + xyz.at(3) + modeTail);
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
    std::cout << " Done" << std::endl;
}