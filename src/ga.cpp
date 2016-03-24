//#include "ga.h"
//#include "nfo2.h"
//
//
////Using genetic algorithm to refine the stable direction of the Mon. matrix
//void gac(gsl_matrix_complex **DAT, gsl_vector_complex* eigenVs, gsl_complex eigenLs)
//{
//
//    cout << "GA" << endl;
//
//    // See if we've been given a seed to use (for testing purposes).  When you
//    // specify a random seed, the evolution will be exactly the same each time
//    // you use that seed number.
//  unsigned int seed = 0;
//
//    // Declare variables for the GA parameters and set them to some default values.
//  int popsize  = 30;
//  int ngen     = 200;
//  float pmut   = 0.01;
//  float pcross = 0.6;
//
//  // Create a phenotype for seven variables.
//  GABin2DecPhenotype map;
//  long double medianV;
//  long double factor = 1e-7;
//
//  //eigenvalue
//  medianV = GSL_REAL(eigenLs);
//  map.add(63, medianV*(1.0-factor), medianV*(1.0+factor));
//  //eigenvector
//  for(int i =0; i<6; i++)
//  {
//    medianV = GSL_REAL(gsl_vector_complex_get(eigenVs, i));
//    map.add(63, medianV*(1.0-factor), medianV*(1.0+factor));
//  }
//
//
//    // Create the template genome using the phenotype map we just made.
//  GABin2DecGenome genome(map, objective);
//  genome.userData(DAT);
//
//    // Now create the GA using the genome and run it.  We'll use sigma truncation
//    // scaling so that we can handle negative objective scores.
//
//  GASimpleGA ga(genome);
//  GASigmaTruncationScaling scaling;
//  ga.populationSize(popsize);
//  ga.nGenerations(ngen);
//  ga.pMutation(pmut);
//  ga.pCrossover(pcross);
//  ga.scaling(scaling);
//  ga.scoreFilename("bog.dat");
//  ga.scoreFrequency(10);
//  ga.flushFrequency(50);
//  ga.evolve(seed);
//
//// Dump the results of the GA to the screen.
//
//  genome = ga.statistics().bestIndividual();
//  cout << "the ga found an optimum at the point (";
//  cout << "ls =  " << genome.phenotype(0) << endl;
//  cout << "vs = (" << genome.phenotype(1) << ", " << genome.phenotype(2) << ", " << genome.phenotype(3) << ", " << genome.phenotype(4) << ", " << genome.phenotype(5) << ", " << genome.phenotype(6) << ")\n\n";
//  cout << "best of generation data are in '" << ga.scoreFilename() << "'\n";
//
//}
//
////Objective function
//float
//objective(GAGenome & c)
//{
//  GABin2DecGenome & genome = (GABin2DecGenome &)c;
//
//  //Eigenvalue
//  gsl_complex eigenLs = gslc_complex(genome.phenotype(0), 0.0);
//  //Eigenvector
//  gsl_vector_complex* eigenVs = gsl_vector_complex_alloc(6);
//  for(int i=0; i<6;i++) gsl_vector_complex_set(eigenVs, i, gslc_complex((double)genome.phenotype(i+1),0.0));
//
//
//  //cout << "eigenLs = " << GSL_REAL(eigenLs) << endl;
//  //cout << "eigenVs = " << endl;
//  //gslc_vector_complex_printf(eigenVs);
//
//
//  gsl_vector_complex* x = gsl_vector_complex_alloc(6);
//  //x = M*eigenVs/eigenLs
//  gslc_matrix_vector_product( (gsl_matrix_complex**) genome.userData(), eigenVs , x, 20);
//  gsl_vector_complex_scale(x, gsl_complex_inverse(eigenLs));
//  //x = x-eigenVs
//  gsl_vector_complex_sub(x, eigenVs);
//
//  //cout << "x = " << endl;
//  //gslc_vector_complex_printf(x);
//
//   //y = ||x||
//  float y = gsl_blas_dznrm2(x);
//
//  cout << "ga. y = " << gsl_blas_dznrm2(x) << ",  " << y << endl;
//
//  //Free memory
//  gsl_vector_complex_free(eigenVs);
//  gsl_vector_complex_free(x);
//
//  return -y;
//}
