/**
 * A fake program that performs fake computations
 * on fake trees and fake sites. It computes superfunctions
 * for each site on one given tree, and does a reduce on the
 * results. Then it iterates again on each tree. 
 * Each rank has a subset of sites.
 */

#include <mpi.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <omp.h>

using namespace std;

int superfunction(unsigned int s, const vector<int> &v)
{
  int res = 0;
  for (auto d: v)
    res += d * s + s % (d+1);
  return res;
}

void fake_openmp_program(int trees, int sites, int threads) 
{
  cout << "Fake computation on " << trees << " trees " <<
      "and " << sites << " sites with " << threads << " threads."<< endl;
  vector<int> bidon_input;
  for (unsigned int i = 0; i < 1000; ++i) {
    bidon_input.push_back(i);
  }
  int globalTreell = 0;
  #pragma omp parallel for num_threads(threads)
  for (unsigned int tree = 0; tree < trees; ++tree) {
    int localTreell = 0;
    for (unsigned int s = 0; s < sites; ++s) {
      localTreell += superfunction(s, bidon_input);      
    }
    #pragma omp critical
    {
      globalTreell += localTreell;
    }
  }
  cout << "result: " <<  globalTreell << endl;
}


void fake_mpi_program(int trees, int sites) 
{

  int rank;
  int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
 
  if (rank == 0) {
    cout << "Fake computation on " << trees << " trees " <<
      "and " << sites << " sites with " << size << " ranks."<< endl;
  }
  
  vector<int> bidon_input;
  for (unsigned int i = 0; i < 1000; ++i) {
    bidon_input.push_back(i);
  }
  int globalTreell = 0;
  int start = (sites / size) * rank;
  int end = min(sites, start + sites / size);
  for (unsigned int tree = 0; tree < trees; ++tree) {
    int localTreell = 0;
    for (unsigned int s = start; s < end; ++s) {
      localTreell += superfunction(s, bidon_input);      
    }
    int temp;
    MPI_Reduce(&localTreell, &temp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    globalTreell += temp;
  }

  if (rank == 0) {
    cout << "result: " <<  globalTreell << endl;
  }

}

bool isMPI() {
  int size = 0; 
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size > 1;
}

int main(int argc, char *argv[])
{
  if (argc != 5) {
    cerr << "Invalid syntax." << endl;
    cerr << "fake_program <trees_number> <sites_number> <parallelization> <openmp_threads>" << endl;
    cerr << "when running with mpi, openmp_threads should be 1" << endl;
    exit(0);
  }

  for (int i = 0; i < argc; ++i) {
    cout << argv[i] << " ";
  }
  cout  << endl;
  int trees = atoi(argv[1]);
  int sites = atoi(argv[2]);
  string implem = argv[3];
  int openmpThreads = atoi(argv[4]);
  int rank = 0;
  auto start = chrono::system_clock::now();
 
  
  if (implem == "mpi") {
    MPI_Init(&argc,&argv);        
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (openmpThreads > 1) {
      cerr << "ERROR: when running with MPI, openmp_threads should be 1" << endl;
      exit(1);
    }
    fake_mpi_program(trees, sites);
    MPI_Finalize();              
  } else if(implem == "openmp") {
    fake_openmp_program(trees, sites, openmpThreads);
  } else {
    cerr << "invalid implem! use openmp or mpi" << endl;
    exit(1);
  }

  if (rank == 0) {
    auto end = chrono::system_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>
                                   (end-start).count();
    cout << "Elapsed time " << elapsed << "ms" << endl; 
  }
  return 0;
}


