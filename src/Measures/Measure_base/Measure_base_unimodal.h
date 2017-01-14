//////////////////////////////////////////////////////////////////////////////////
//    Copyright (C) 2016,  Constantinos Tsirogiannis.  Email: tsirogiannis.c@gmail.com
//
//    This file is part of PhyloMeasures.
//
//    PhyloMeasures is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    PhyloMeasures is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with PhyloMeasures.  If not, see <http://www.gnu.org/licenses/>
//////////////////////////////////////////////////////////////////////////////////

#ifndef MEASURE_BASE_UNIMODAL_H
#define MEASURE_BASE_UNIMODAL_H

#include<vector>
#include<cmath>

namespace PhylogeneticMeasures {

template< class KernelType >
struct Measure_base_unimodal
{
  typedef KernelType                                        Kernel;
  typedef Measure_base_unimodal<Kernel>                     Self;
  typedef typename Kernel::Number_type                      Number_type;
  typedef typename Kernel::Distribution_type                Distribution_type;
  typedef typename Kernel::Exception_type                   Exception_type;
  typedef typename Kernel::Exception_functor                Exception_functor;
  typedef typename Kernel::Uniform_sampler                  Uniform_sampler;  
  typedef typename Kernel::Sequential_sampler               Sequential_sampler;  
  typedef typename Kernel::Incremental_Monte_Carlo_handler  Incremental_Monte_Carlo_handler;  

 public: 
 
  Measure_base_unimodal():_distribution(Kernel::UNIFORM_FIXED_SIZE),_seed(0){}

  void set_probability_distribution( Distribution_type distrib )
  { _distribution = distrib; }

  Distribution_type probability_distribution() const
  { return _distribution; }

  Self& operator=(const Self& d)
  {
    _distribution = d.probability_distribution();
    return *this;
  }

  int seed() const
  { return _seed;}

  void set_seed(int seed)
  { _seed = seed; }
  
 protected:

  // Function that reads a list of integers from a file,
  // assumed to be samples sizes that are used as input
  // for computing the expectation and deviation of a measure
  // on a given tree. 
  template < class TreeType, class OutputIterator >
  void _read_sample_sizes_from_file(char *filename, TreeType &tree, OutputIterator ot);
  
  // Input:  A tree and a range of iterators that indicate a list of species names (in std::string format).
  // Output: The value of the current measure for this set of species.
  template <class TreeType, class RangeIterator, class Measure>    
  Number_type _list_query(TreeType &tree, RangeIterator rbegin, RangeIterator rend, Measure &msr);
  
  // Input: A tree and a txt file that stores a list of species names, which constitute a subset 
  // of the species (leaf nodes) in the tree and which appear in random order.
  // Output: The value of the current measure for this set of species.
  template <class TreeType, class Measure>
  Number_type _list_query(TreeType &tree, char* filename, Measure &msr);


  void read_csv_matrix(const char *filename, std::vector<std::string> &names, 
                       std::vector< std::vector<bool> > &matrix);

  // Input: A tree, an array with the species names in the tree and a matrix such that: each column corresponds 
  // to one of these species and each row indicates a sample of these species for we want to compute the
  // distance measure. A certain species is considered as part of the i-th sample if in the i-th row 
  // of the matrix there is a '1' at the column that corresponds to this species (otherwise a '0').

  // The fifth argument is a boolean value that indicates if we just want to compute the value of 
  // the measure for the indicated sample of species, or if we want to "standardise" this value;
  // that is to subtract from this value the mean of the measure for all samples of the same size
  // and then divide by the standard deviation.

  // The sixth argument is an output iterator of the type std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the (standardised) value of the measure for the sample that is described in the i-th row of the input matrix.

  template < class TreeType, class Measure, class OutputIterator>
  int _matrix_query( TreeType &tree, std::vector<std::string> &names, 
                     std::vector< std::vector<bool> > &matrix, 
                     Measure &msr, bool standardised, OutputIterator ot );
  
  // Input: A csv file which stores a matrix where each column corresponds to a species of the tree
  // and each row indicates a sample of these species for which we want to compute the
  // distance measure. A certain species is considered as part of the i-th sample if in the i-th row 
  // of the matrix there is a '1' at the column that corresponds to this species (otherwise there is a '0').

  // The fourth argument is a boolean value that indicates if we just want to compute the value of 
  // the measure for the indicated sample of species, or if we want to "standardise" this value;
  // that is to subtract from this value the mean of the measure for all samples of the same size
  // and then divide by the standard deviation.

  // The fifth argument is an output iterator of the type std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the (standardised) value of the measure for the sample that is described in the i-th row of the input matrix.
  template <class TreeType, class Measure, class OutputIterator>
  int _csv_matrix_query( TreeType &tree, char *filename, Measure &msr, bool standardised, OutputIterator ot );


  template < class TreeType, class Measure, class OutputIterator>
  int _matrix_query_Poisson_binomial_fixed_size( TreeType &tree, std::vector<std::string> &names, 
                                                 std::vector< std::vector<bool> > &matrix,
                                                 Measure &msr, bool standardised, OutputIterator ot );

  template < class TreeType, class Measure, class OutputIterator>
  int _csv_matrix_query_Poisson_binomial_fixed_size( TreeType &tree, char *filename, Measure &msr, 
                                                      bool standardised, OutputIterator ot );

  // Uses Monte-Carlo methods
  template < class TreeType, class Measure, class OutputIterator>
  int _matrix_query_sequential_fixed_size( TreeType &tree, std::vector<std::string> &names, 
                                           std::vector< std::vector<bool> > &matrix,
                                           Measure &msr, bool standardised, 
                                           OutputIterator ot, int repetitions=1000);

  // Uses Monte-Carlo methods
  template < class TreeType, class Measure, class OutputIterator>
  int _csv_matrix_query_sequential_fixed_size( TreeType &tree, char *filename, Measure &msr, 
                                               bool standardised, OutputIterator ot, int repetitions=1000);
  
  template < class TreeType, class Measure, class OutputIterator>
  int _matrix_query_weighted( TreeType &tree, std::vector<std::string> &names, 
                               std::vector< std::vector<bool> > &matrix,
                               Measure &msr, bool standardised, OutputIterator ot, int repetitions=1000 )
  {
    if(msr.probability_distribution() == Kernel::POISSON_BINOMIAL_FIXED_SIZE)
      return _matrix_query_Poisson_binomial_fixed_size( tree, names, matrix, msr, standardised, ot);
    else if(msr.probability_distribution() == Kernel::SEQUENTIAL_FIXED_SIZE)
      return _matrix_query_sequential_fixed_size( tree, names, matrix, msr, standardised, ot, repetitions);
    else
      return 0;
  }

  template < class TreeType, class Measure, class OutputIterator>
  int _csv_matrix_query_weighted( TreeType &tree, char *filename, Measure &msr, 
                                   bool standardised, OutputIterator ot, int repetitions=1000 )
  {
    if(msr.probability_distribution() == Kernel::POISSON_BINOMIAL_FIXED_SIZE)
      return _csv_matrix_query_Poisson_binomial_fixed_size( tree, filename, msr, standardised, ot);
    else if(msr.probability_distribution() == Kernel::SEQUENTIAL_FIXED_SIZE)
      return _csv_matrix_query_sequential_fixed_size( tree, filename, msr, standardised, ot, repetitions);
    else
      return 0;
  }

  /////////////////////////////////////
  // Functions that compute p-values // 
  /////////////////////////////////////

  // Uses Monte-Carlo methods
  template < class TreeType, class Measure, class OutputIterator>
  int _pvalues_query_uniform_fixed_size( TreeType &tree, std::vector<std::string> &names, 
                                         std::vector< std::vector<bool> > &matrix,
                                         Measure &msr, OutputIterator ot, int repetitions=1000);

  // Uses Monte-Carlo methods
  template < class TreeType, class Measure, class OutputIterator>
  int _csv_pvalues_query_uniform_fixed_size( TreeType &tree, const char *filename, Measure &msr, 
                                             OutputIterator ot, int repetitions=1000);

  // Uses Monte-Carlo methods
  template < class TreeType, class Measure, class OutputIterator>
  int _pvalues_query_sequential_fixed_size( TreeType &tree, std::vector<std::string> &names, 
                                            std::vector< std::vector<bool> > &matrix,
                                            Measure &msr, OutputIterator ot, int repetitions=1000);

  // Uses Monte-Carlo methods
  template < class TreeType, class Measure, class OutputIterator>
  int _csv_pvalues_query_sequential_fixed_size( TreeType &tree, const char *filename, Measure &msr, 
                                                OutputIterator ot, int repetitions=1000);

  ////////////////////////////////////
  // Functions that compute moments // 
  ////////////////////////////////////

  // Uses Monte-Carlo methods
  template < class Measure, class OutputIterator>
  void _compute_moments_sequential_fixed_size( Measure &msr, int sample_size, 
                                               OutputIterator ot_a, OutputIterator ot_b, 
                                               int repetitions=1000);

 private:
    
  Distribution_type        _distribution;
  int                     _seed;

}; // Measure_base_unimodal

} // namespace PhylogeneticMeasures

#include"Measure_base_unimodal_impl.h"

#endif //MEASURE_BASE_UNIMODAL_H
