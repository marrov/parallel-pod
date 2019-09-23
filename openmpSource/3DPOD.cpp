#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <sys/time.h>
#include <omp.h>

#define MSIZE 381600 // Rows of matrix (number of points)
#define TSIZE 5      // Size of time data (number of snapshots)
#define VSIZE 3      // Size of the variable (1 if scalar, 3 if vector)
#define NSIZE 5      // Size of output POD modes (number of modes to write)

using namespace Eigen;

int main()
{
    float runTime=0;
    struct timeval start;
    struct timeval end;

    MatrixXd m = MatrixXd::Zero(MSIZE * VSIZE, TSIZE);
    MatrixXd pm = MatrixXd::Zero(TSIZE, TSIZE);

    // READING INPUT FILES

    gettimeofday(&start,NULL);
    omp_set_num_threads(4);
#pragma omp parallel
#pragma omp for
    for (size_t k = 0; k < TSIZE; k++)
    {
        std::string dir = "input/U/";
        std::string fname = "U_" + std::to_string(k + 1) + ".dat";
        std::ifstream file(dir + fname);

        if (file.is_open())
        {
            for (size_t i = 0; i < MSIZE; i++)
            {
                for (size_t j = 0; j < VSIZE; j++)
                {
                    file >> m(i + MSIZE * j, k);
                }
            }
//            std::cout << "Finished reading " << fname << std::endl;
            std::cout << "Finished reading " << fname <<" by thread" << omp_get_thread_num() << std::endl;
            file.close();
        }
        else
        {
            std::cout << " Unable to open file" << std::endl;
        }
    }
    gettimeofday(&end,NULL);
    runTime=(end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec);
    std::cout << "Reading files takes "<< runTime/1000000 << "s"<< std::endl;


    // COMPUTING NORMALISED PROJECTION MATRIX

    gettimeofday(&start,NULL);
#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < TSIZE; i++)
    {
        for (size_t j = 0; j < TSIZE; j++)
        {
            pm(i, j) = (1.0 / TSIZE) * (m.col(i).dot(m.col(j).transpose()));
        }
    }
    gettimeofday(&end,NULL);
    runTime=(end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec);
    std::cout << "computing projection matrix takes "<< runTime/1000000 << "s"<< std::endl;

    /*
    std::cout << "Normalised rojection Matrix:" << std::endl
              << std::setprecision(12) << pm << std::endl;
    */

    // COMPUTING SORTED EIGENVALUES AND EIGENVECTORS

    gettimeofday(&start,NULL);
    SelfAdjointEigenSolver<MatrixXd> eigensolver(pm);
    if (eigensolver.info() != Success)
        abort();

    VectorXd eigval = eigensolver.eigenvalues().reverse();
    MatrixXd eigvec = eigensolver.eigenvectors().rowwise().reverse();
    gettimeofday(&end,NULL);
    runTime=(end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec);
    std::cout << "eigen-decomposition takes "<< runTime/1000000 << "s"<< std::endl;

    /*
    std::cout << "Sorted eigenvalues of Projection Matrix:\n"
              << eigval << std::endl;
    std::cout << "Matrix whose columns are the sorted eigenvectors of Projection Matrix:\n"
              << eigvec << std::endl;
    */

    // WRITING SORTED EIGENVALUES

    std::ofstream writeEigval("chronos/A.txt");
    if (writeEigval.is_open())
    {
        writeEigval << std::scientific << std::setprecision(10) << eigval;
        writeEigval.close();
    }

    // COMPUTING POD MODES AND WRITING CHRONOS FOR EACH MODE

    MatrixXd podx = MatrixXd::Zero(MSIZE, TSIZE);
    MatrixXd pody = MatrixXd::Zero(MSIZE, TSIZE);
    MatrixXd podz = MatrixXd::Zero(MSIZE, TSIZE);

    gettimeofday(&start,NULL);
#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < TSIZE; i++)
    {
        std::string chronos = "chronos/chronos_" + std::to_string(i) + ".dat";
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
    gettimeofday(&end,NULL); 
    runTime=(end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec);
    std::cout << "computing POD modes takes "<< runTime/1000000 << "s"<< std::endl;

    // WRITING POD MODES
    gettimeofday(&start,NULL);
#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < NSIZE; i++)
    {
        std::string modex = "mode/mode_Ux_" + std::to_string(i) + ".dat";
        std::string modey = "mode/mode_Uy_" + std::to_string(i) + ".dat";
        std::string modez = "mode/mode_Uz_" + std::to_string(i) + ".dat";

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
    gettimeofday(&end,NULL);
    runTime=(end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec);
    std::cout << "writing POD modes takes "<< runTime/1000000 << "s"<< std::endl;
}
