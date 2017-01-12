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

#ifndef CORE_ANCESTOR_COST_H
#define CORE_ANCESTOR_COST_H

#include<string>
#include<vector>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<set>
#include<queue>
#include<map>
#include<cmath>


namespace PhylogeneticMeasures {

template< class KernelType >
struct Core_ancestor_cost: public KernelType::Measure_base_unimodal
{
  typedef KernelType                                       Kernel;
  typedef typename Kernel::Measure_base_unimodal          Base;
  typedef typename Kernel::Number_type                    Number_type;
  typedef typename Kernel::Unimodal_tree                  Tree_type;
  typedef typename Tree_type::Node_type                   Node_type;
  typedef typename Kernel::Numeric_traits                 Numeric_traits;
  typedef typename Numeric_traits::Protected_number_type  Protected_number_type;
  typedef typename Numeric_traits::Square_root            Square_root;
  typedef typename Numeric_traits::Ceiling                Ceiling;
  typedef typename Numeric_traits::To_double              To_double;

  typedef typename Kernel::Exception_type                 Exception_type;
  typedef typename Kernel::Exception_functor              Exception_functor; 

  struct Data_type
  {
    Number_type chi;
  };

  public:

  Core_ancestor_cost(Tree_type &tree, Number_type chi=0.501)
  { 
    p_tree = &tree;

    if( chi <= Number_type(0.5) || chi > Number_type(1.0) )
    {
      std::string msg;
      msg += " Invalid value of parameter chi. The value of chi must belong to the interval (0.5,1.0] .\n";
      Exception_type excp;
      excp.get_error_message(msg);
      Exception_functor excf;
      excf(excp);
    }

    _chi = chi;
  }

  Tree_type& tree(void)
  { return *p_tree;}

  Tree_type* tree_pointer(void)
  { return p_tree;}

  Data_type auxiliary_data() const
  {
    Data_type dt;
    dt.chi = _chi; 
    return dt;
  }

  Number_type reference_value()
  {
    Number_type val(-1.0);

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      if( p_tree->node(i).distance > Number_type(0.0) && 
          ( val <= Number_type(0.0) || p_tree->node(i).distance < val ) )
        val = p_tree->node(i).distance;

    return val;
  }

  void set_auxiliary_data(const Data_type &dt)
  {
    if( dt.chi <= Number_type(0.5) || dt.chi > Number_type(1.0) )
    {
      std::string msg;
      msg += " Invalid value of parameter chi. The value of chi must belong to the interval (0.5,1.0] .\n";
      Exception_type excp;
      excp.get_error_message(msg);
      Exception_functor excf;
      excf(excp);
    }

    _chi = dt.chi;

  } // set_auxiliary_data(const Data_type &dt)


  Number_type chi()
  { return _chi;}

  template< class RangeIterator >
  Number_type operator()( RangeIterator rbegin, RangeIterator rend,
                          int min_index, int max_index );


  template < class OutputIterator >
  void incremental_operator( std::vector<int> &sample,
                              std::vector<int> &sample_sizes, OutputIterator ot );

  Number_type update_marked_subtree(int new_node_index, int rchi);


  template< class RangeIterator >
  Number_type slow_operator( RangeIterator rbegin, RangeIterator rend,
                             int min_index, int max_index);
			  
						    
  // Input:  A range of iterators that indicate a list of species names (in std::string format).
  // Output: The value of the current measure for this set of species.
  template <class RangeIterator>    
  Number_type list_query(RangeIterator rbegin, RangeIterator rend)
  {return this->_list_query(*p_tree, rbegin, rend, *this);}
  
  // Input: A txt file that stores a list of species names, which constitute a subset 
  // of the species (leaf nodes) in the tree and which appear in random order.
  // Output: The value of the current measure for this set of species.
  Number_type list_query(char* filename)
  {return this->_list_query(*p_tree, filename, *this);}


  // Input: A vector with the species names from the tree, and a matrix such that: each column 
  // corresponds to one of these species, and each row indicates a sample of these species 
  // for which we want to compute the distance measure. A certain species is considered as 
  // part of the i-th sample if in the i-th row of the matrix there is a '1' at the column
  // that corresponds to this species (otherwise there is a '0').

  // The last argument is an output iterator of the type std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the value of the measure for the sample that is described in the i-th row of the input matrix.
  template <class OutputIterator>
  int matrix_query_basic( std::vector<std::string> &names, 
                          std::vector< std::vector<bool> > &matrix, OutputIterator ot )
  {return this->_matrix_query(*p_tree, names, matrix, *this, false, ot);}


  // Same as the function above except that the returned values distance values have been "standardised":
  // from each value we have subtracted the mean value of the measure among all samples of leaves of the
  // same size, and we have divided the result by the deviation.
  template <class OutputIterator>
  int matrix_query_standardised(std::vector<std::string> &names, 
                                std::vector< std::vector<bool> > &matrix, OutputIterator ot)
  {return this->_matrix_query(*p_tree, names, matrix, *this, true, ot);}



  template <class OutputIterator>
  int matrix_query_weighted_basic(std::vector<std::string> &names, 
                                  std::vector< std::vector<bool> > &matrix, 
                                  OutputIterator ot, int repetitions=1000 )
  {
    if(this->probability_distribution() != Kernel::SEQUENTIAL_FIXED_SIZE)
    {
      std::string exception_msg;
      exception_msg += " The distribution of the CAC object should be set to"; 
      exception_msg += " Kernel::SEQUENTIAL_FIXED_SIZE .";
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    } 

    return this->_matrix_query_weighted(*p_tree, names, matrix, *this, false, ot, repetitions);
  }


  template <class OutputIterator>
  int matrix_query_weighted_standardised(std::vector<std::string> &names, 
                                         std::vector< std::vector<bool> > &matrix, 
                                         OutputIterator ot, int repetitions=1000 )
  {
    if(this->probability_distribution() != Kernel::SEQUENTIAL_FIXED_SIZE)
    {
      std::string exception_msg;
      exception_msg += " The distribution of the CAC object should be set to"; 
      exception_msg += " Kernel::SEQUENTIAL_FIXED_SIZE .";
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    } 

    return this->_matrix_query_weighted(*p_tree, names, matrix, *this, true, ot, repetitions);
  }
  
  // Input: A csv file which stores a matrix where each column corresponds to a species of the tree
  // and each row indicates a sample of these species for which we want to compute the
  // distance measure. A certain species is considered as part of the i-th sample if in the i-th row 
  // of the matrix there is a '1' at the column that corresponds to this species (otherwise there is a '0').

  // The second argument is an output iterator of the type std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the value of the measure for the sample that is described in the i-th row of the input matrix.
  template <class OutputIterator>
  int csv_matrix_query_basic(char *filename, OutputIterator ot )
  {return this->_csv_matrix_query(*p_tree, filename, *this, false, ot);}


  // Same as the function above except that the returned values distance values have been "standardised":
  // from each value we have subtracted the mean value of the measure among all samples of leaves of the
  // same size, and we have divided the result by the deviation.
  template <class OutputIterator>
  int csv_matrix_query_standardised(char *filename, OutputIterator ot)
  {return this->_csv_matrix_query(*p_tree, filename, *this, true, ot);}


  template <class OutputIterator>
  int csv_matrix_query_weighted_standardised(char *filename, OutputIterator ot, int repetitions=1000 )
  {
    if(this->probability_distribution() != Kernel::SEQUENTIAL_FIXED_SIZE)
    {
      std::string exception_msg;
      exception_msg += " The distribution of the CAC object should be set to"; 
      exception_msg += " Kernel::SEQUENTIAL_FIXED_SIZE .";
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    } 

    return this->_csv_matrix_query_weighted(*p_tree, filename, *this, true, ot, repetitions);
  }

  template <class OutputIterator>
  int csv_matrix_query_weighted_basic(char *filename, OutputIterator ot )
  {return this->_csv_matrix_query(*p_tree, filename, *this, false, ot);}

  /////////////////////////////////////
  // Functions that compute p-values //
  /////////////////////////////////////

  // Uses Monte-Carlo methods
  template <class OutputIterator>
  int pvalues_query_uniform_fixed_size( std::vector<std::string> &names, 
                                        std::vector< std::vector<bool> > &matrix,
                                        OutputIterator ot, int repetitions=1000)
  { return this->_pvalues_query_uniform_fixed_size( *p_tree, names, matrix,*this, ot, repetitions); }

  // Uses Monte-Carlo methods
  template <class OutputIterator>
  int csv_pvalues_query_uniform_fixed_size(const char *filename, OutputIterator ot, int repetitions=1000)
  { return this->_csv_pvalues_query_uniform_fixed_size( *p_tree, filename, *this, ot, repetitions); }

  // Uses Monte-Carlo methods
  template <class OutputIterator>
  int pvalues_query_sequential_fixed_size( std::vector<std::string> &names, 
                                           std::vector< std::vector<bool> > &matrix,
                                           OutputIterator ot, int repetitions=1000)
  { return this->_pvalues_query_sequential_fixed_size(*p_tree, names, matrix, *this, ot, repetitions);}

  // Uses Monte-Carlo methods
  template <class OutputIterator>
  int csv_pvalues_query_sequential_fixed_size(const char *filename, OutputIterator ot, int repetitions=1000)
  { return this->_csv_pvalues_query_sequential_fixed_size( *p_tree, filename, *this, ot, repetitions); }


  // Function that reads a list of integers from a file,
  // assumed to be samples sizes that are used as input
  // for computing the expectation and deviation of a measure
  // on a given tree. 
  template < class OutputIterator >
  void read_sample_sizes_from_file(char *filename, OutputIterator ot)
  { this->_read_sample_sizes_from_file(filename, *p_tree, ot);}


  // Computes all together the probability values f(x) = \binom{x}{r}/\binom{s}{r}

  void compute_all_hypergeometric_probabilities_a( int sample_size, int number_of_leaves);

  void compute_all_hypergeometric_probabilities_b( int sample_size, int number_of_leaves);

  Protected_number_type hypergeom_a( int x )
  {
    if( x < _sample_size || x > _number_of_leaves )
      return Protected_number_type(Number_type(0.0));

    if( x == _number_of_leaves )
      return Protected_number_type(Number_type(1.0));

    return _hypergeom_a[x-_sample_size];
  }

  Protected_number_type hypergeom_b( int x )
  {
    if( x < _number_of_leaves - _sample_size || x > _number_of_leaves )
      return Protected_number_type(Number_type(0.0));

    if( x == _number_of_leaves )
      return Protected_number_type(Number_type(1.0));

    return _hypergeom_b[x-_number_of_leaves+_sample_size];
  }

  template<class OutputIterator>
  void compute_k_binomial_coefficients(int k, OutputIterator ot)
  {
    Protected_number_type coeff(1.0);

    *ot++ = coeff;

    for(int i=0; i<k; i++)
    {
      coeff = coeff*Protected_number_type(Number_type(k-i))/
                    Protected_number_type(Number_type(i+1));

      *ot++ = coeff;
    }
  }

  Protected_number_type compute_node_probability
  (int number_of_subtree_leaves, int sample_size, bool silent = true);

  template<class OutputIterator>
  void compute_all_root_path_costs(OutputIterator ot);

  template<class OutputIterator>
  void compute_first_k_raw_moments(int k, int sample_size, OutputIterator ot);


  template<class OutputIterator>
  void compute_first_k_raw_moments_protected(int k, int sample_size, OutputIterator ot);

  template <class OutputIterator>
  void compute_first_k_centralised_moments( int k, int sample_size, OutputIterator ot);


  Number_type compute_expectation( int sample_size )
  {
    if( sample_size < 0 || sample_size > p_tree->number_of_leaves() )
    {
      std::string exception_msg;
      exception_msg += " Request to compute expectation with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }


    if(this->probability_distribution() == Kernel::UNIFORM_FIXED_SIZE)
    {
      std::vector<Number_type> mean;
      compute_first_k_raw_moments(1, sample_size, std::back_inserter(mean));
      return mean[0];
    }
    else if(this->probability_distribution() == Kernel::SEQUENTIAL_FIXED_SIZE)
    {
      if(sample_size > _Sequential_exps.size()-1 || _Sequential_exps.size() == 0)
      {
        _Sequential_exps.clear();
        _Sequential_devs.clear();
       
        this->_compute_moments_sequential_fixed_size(*this, sample_size, 
                                                      std::back_inserter(_Sequential_exps),
                                                      std::back_inserter(_Sequential_devs),1000);
      }

      return _Sequential_exps[sample_size];

    } // else if(this->probability_distribution() == Kernel::SEQUENTIAL_FIXED_SIZE)
    else
      return Number_type(-1.0);

  } // compute_expectation(...)

  Number_type compute_variance( int sample_size )
  {
    if( sample_size < 0 || sample_size > p_tree->number_of_leaves() )
    {
      std::string exception_msg;
      exception_msg += " Request to compute variance with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }


    if(this->probability_distribution() == Kernel::UNIFORM_FIXED_SIZE)
    {
      std::vector<Number_type> moments;

      compute_first_k_centralised_moments(2, sample_size, std::back_inserter(moments));
      return moments[1];
    }
    else if(this->probability_distribution() == Kernel::SEQUENTIAL_FIXED_SIZE)
    {
      if(sample_size > _Sequential_exps.size()-1 || _Sequential_exps.size() == 0)
      {
        _Sequential_exps.clear();
        _Sequential_devs.clear();
       
        this->_compute_moments_sequential_fixed_size(*this, sample_size, 
                                                      std::back_inserter(_Sequential_exps),
                                                      std::back_inserter(_Sequential_devs),1000);
      }

      return _Sequential_devs[sample_size]*_Sequential_devs[sample_size];

    } // else if(this->probability_distribution() == Kernel::SEQUENTIAL_FIXED_SIZE)
    else
      return Number_type(-1.0);

  } // compute_variance(...)


  Number_type compute_deviation( int sample_size )
  { 
    if( sample_size < 0 || sample_size > p_tree->number_of_leaves() )
    {
      std::string exception_msg;
      exception_msg += " Request to compute deviation with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(this->probability_distribution() == Kernel::UNIFORM_FIXED_SIZE)
    {
      Number_type variance = compute_variance(sample_size);

      if( variance < Number_type(0.0) ) 
        return Number_type(0.0);

      return Square_root()(variance);
    }
    else if(this->probability_distribution() == Kernel::SEQUENTIAL_FIXED_SIZE)
    {
      if(sample_size > _Sequential_exps.size()-1 || _Sequential_exps.size() == 0)
      {
        _Sequential_exps.clear();
        _Sequential_devs.clear();
       
        this->_compute_moments_sequential_fixed_size(*this, sample_size, 
                                                      std::back_inserter(_Sequential_exps),
                                                      std::back_inserter(_Sequential_devs),1000);
      }

      return _Sequential_devs[sample_size];

    } // else if(this->probability_distribution() == Kernel::SEQUENTIAL_FIXED_SIZE)
    else
      return Number_type(-1.0);


  } // compute_deviation(...)

  private:

  Tree_type                  *p_tree; // Stores a pointer to a phylogenetic tree object.
  std::vector<Protected_number_type>   _hypergeom_a, _hypergeom_b; 
                                      // Store values of two hypergeometric functions
  Number_type                _chi;
  std::vector<Number_type>   _Sequential_exps, _Sequential_devs;
  int                        _sample_size;
  int                        _number_of_leaves;

}; // struct Core_ancestor_cost

} // namespace PhylogeneneticMeasures

#include "Core_ancestor_cost_impl.h"

#endif // CORE_ANCESTOR_COST_H
