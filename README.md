# parallel-pod
PDC Summer School project to port a python POD code into C and make it parallel

# Introduction:
The evolution of computational power, and the increasing availability of HPC resources to researchers, allows for pushing the envelope regarding scientific data generation and the subsequent data processing. In the area of fluid dynamics the Proper Orthogonal Decomposition (POD) method is increasingly being used to isolate constituting behaviors in fluid flows in the post-processing stage. The POD method is applicable to both numerically generated as well as experimentally generated datasets. Considering contemporary simulations and numerical experiments, the generated data set can be large. Due to the presence of large datasets of interest for post-processing effective and efficient algorithms and their implementations are of high interest.
 
# Project description:
At present the POD method implementation in use in our close proximity is in Python and of sequential type. Typically, applying this implementation to commonly encountered datasets require around 24 wall-clock hours to complete. It is believed that a performance-improvement can be achieved by, first of all, porting the implementation to a compiled language (such as, C/C++ or Fortran). The current implementation consists of several nested for-loops which should allow themselves for efficient parallelization. Secondly, as the POD method is in practice a Singular Value Decomposition (SVD), it is believed that both an own implementation as well as exploration of available linear algebra packages can be included in the project.
 
Furthermore, the method used to generate input data to the POD algorithm is currently generating several small files on persistent storage. It is not expected to reduce the execution time to a significant extent, but the data loading (I/O) could also be considered for parallelization in this project. If time allows.
 
# Project members:
Marc Rovira (marrs@kth.se), Yazhou Shen (yashen@kth.se), Kristian Rönnberg (kriron@kth.se)
 
# Intended outcome:
We are currently all working on CFD at the Mechanics department at KTH. All of as either are or are intending to apply POD in our research. POD is applied to datasets generated with OpenFOAM. We currently have good datasets, and a reference implementation, which we can be used as starting points for the project. Kristian is not located at KTH more than once per week, so we intend to exercise version control and usage of GitHub to store and manage our progress. Furthermore, our intention is first to port the implementation from Python to C/C++, which will allow for a benchmarking of a compiled vs. interpreted language. The second step is to parallelize the code using MPI. The reason for choosing MPI is that we are currently working extensively with OpenFOAM, which is parallelized using MPI, and thus consider this to be a good chance to get a glimpse of how this works. Hence, SMP is not considered at the moment. As a third, and optional step, an OpenMP or CUDA-based implementation could be considered.
 
In addition to learning more about parallel programming and MPI, a working implementation that will facilitate our future research is the desired outcome. Experience using git and collaborating on a publicly shared and off-site hosted source code repository, such as GitHub, is also of relevance.

# To-do list:
- [ ] Make list of how the input data should be read (format, etc) and how it flows along the code
- [X] (M) Make verification/test case (for RAW input) with run script
- [ ] (M) Profile the python verification case 
- [ ] Understand the legacy python code and write the pseudo code version
- [ ] WIP decide a coding style for C++
- [ ] OPT define a validation case
- [ ] Identify linear algebra packages in C++ (serial and parallel)

Key: OPT: optional, WIP: work in progess, (M): task for Marc
