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

#ifndef MEAN_NEAREST_TAXON_DISTANCE_H
#define MEAN_NEAREST_TAXON_DISTANCE_H

#include<vector>
#include<map>

namespace PhylogeneticMeasures {

template< class KernelType >
struct Mean_nearest_taxon_distance: public KernelType::Measure_base_unimodal
{
  typedef KernelType                                           Kernel;
  typedef typename Kernel::Measure_base_unimodal               Base;
  typedef typename Kernel::Number_type                         Number_type;
  typedef typename Kernel::Mean_nearest_taxon_distance_tree    Tree_type;
  typedef typename Tree_type::Node_type                        Node_type;
  typedef typename Kernel::Numeric_traits                      Numeric_traits;
  typedef typename Numeric_traits::Square_root                 Square_root;
  typedef typename Kernel::Distribution_type                   Distribution_type;
  typedef typename Kernel::Edge_relation_type                  Edge_relation_type;

  typedef typename Kernel::Poisson_binomial_moments_Mean_nearest_taxon_distance 
                                                                 Poisson_moments_functor;

  typedef typename Kernel::Exception_type                      Exception_type;
  typedef typename Kernel::Exception_functor                   Exception_functor;

  struct Data_type{};

 public:

  Mean_nearest_taxon_distance(Tree_type &tree)
  { 
    p_tree = &tree;
    _max_subtree_path_costs.assign(p_tree->number_of_nodes(),Number_type(0.0));
  }

  Tree_type& tree(void)
  { return *p_tree;}

  Tree_type* tree_pointer(void)
  { return p_tree;}

  Data_type auxiliary_data() const
  { return Data_type();}

  void set_auxiliary_data(const Data_type &dt)
  {}

 private:

  Number_type _compute_subtree_min_values( Tree_type &tree, int current_index );

  void _compute_rest_tree_min_values( Tree_type &tree, int current_index );

  Number_type two_edge_pr( int se, int sl , Edge_relation_type er )
  {
    switch(er)
    {
      case Kernel::OFFSPRING:    return hypergeom_minus_one(_number_of_leaves-se);
      case Kernel::ANCESTOR:     return hypergeom_minus_one(_number_of_leaves-sl);
      case Kernel::INDEPENDENT:  return hypergeom_minus_two(_number_of_leaves-sl-se);
    }

    return Number_type(-1.0);
	
  } //two_edge_pr(int se, int sl , Edge_relation_type er )


  template< class OutputIterator >
  void _compute_subtree_sums( int index, Number_type& sum_of_products, OutputIterator ot,
                              Number_type &sum_subtree, Number_type &sum_subtract );


  void _compute_subtree_sums(Number_type &sum_subtree, Number_type &sum_subtract)
  {
    Node_type root = p_tree->root();

    for( int i = 0; i < root.children.size(); i++ )
    {
      std::vector< std::pair<Number_type, int> > subtree_leaves;
      Number_type sum_of_products(0.0);

      _compute_subtree_sums( root.children[i], sum_of_products,
                             std::back_inserter(subtree_leaves),
                             sum_subtree, sum_subtract );

      subtree_leaves.clear();
    }

  } // void _compute_subtree_sums( ... )

 public:


  Number_type reference_value()
  {
    Number_type val(-1.0);

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      if( p_tree->node(i).distance > Number_type(0.0) && 
          ( val <= Number_type(0.0) || p_tree->node(i).distance < val ) )
        val = p_tree->node(i).distance;

    return val;
  }

  template< class RangeIterator >
  Number_type operator()( RangeIterator rbegin, RangeIterator rend,
                          int min_index, int max_index );


  template< class OutputIterator >
  void incremental_operator_non_ultrametric
  ( std::vector<int> &sample, std::vector<int> &sample_sizes, OutputIterator ot );

  template < class OutputIterator >
  void incremental_operator_ultrametric
  ( std::vector<int> &sample, std::vector<int> &sample_sizes, OutputIterator ot );

  template < class OutputIterator >
  void incremental_operator( std::vector<int> &sample,
                             std::vector<int> &sample_sizes, OutputIterator ot )
  {
    if(p_tree->is_ultrametric())
      incremental_operator_ultrametric(sample,sample_sizes,ot);
    else
      incremental_operator_non_ultrametric(sample,sample_sizes,ot);
  }

  Number_type update_total_cost_ultrametric(int &intersection_index, int new_node_index);

  void initialize_max_subtree_path_costs(int index);

  void update_shortest_path_costs(int &intersection_index,int new_node_index, Number_type &total_dist);

  template< class OutputIterator >  
  void find_new_nearest_neighbours(Number_type dist, int index, OutputIterator ot, 
                                   Number_type &distance_difference);

  void update_max_subtree_path_costs(int index);
						    
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

  // Function that reads a list of integers from a file,
  // assumed to be samples sizes that are used as input
  // for computing the expectation and deviation of a measure
  // on a given tree. 
  template < class OutputIterator >
  void read_sample_sizes_from_file(char *filename, OutputIterator ot)
  { this->_read_sample_sizes_from_file(filename, *p_tree, ot);}

  ////////////////////////////////////////////////////
  // Query functions for non-weighted computations. //
  ////////////////////////////////////////////////////

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
  { return this->_matrix_query(*p_tree, names, matrix, *this, false, ot);}


  // Same as the function above except that the returned values distance values have been "standardised":
  // from each value we have subtracted the mean value of the measure among all samples of leaves of the
  // same size, and we have divided the result by the deviation.
  template <class OutputIterator>
  int matrix_query_standardised(std::vector<std::string> &names, 
                                 std::vector< std::vector<bool> > &matrix, OutputIterator ot )
  {
    if(p_tree->is_ultrametric() == true)
      return this->_matrix_query(*p_tree, names, matrix, *this, true, ot);
    else
      return -1;
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
  { return this->_csv_matrix_query(*p_tree, filename, *this, false, ot);}


  // Same as the function above except that the returned values distance values have been "standardised":
  // from each value we have subtracted the mean value of the measure among all samples of leaves of the
  // same size, and we have divided the result by the deviation.
  template <class OutputIterator>
  int csv_matrix_query_standardised(char *filename, OutputIterator ot )
  {
    if(p_tree->is_ultrametric() == true)
      return this->_csv_matrix_query(*p_tree, filename, *this, true, ot);
    else
      return -1;
  }

  //////////////////////////////////////////////////////////
  // Query functions for weighted computations.           //
  // The standardised versions of these measures          //
  // are calculated based on the probability              //
  // distribution which is stored as a flag               //
  // in this measure object (probability_distribution()). // 
  //////////////////////////////////////////////////////////

  // Input: A vector with the species names from the tree, and a matrix such that: each column 
  // corresponds to one of these species, and each row indicates a sample of these species 
  // for which we want to compute the distance measure. A certain species is considered as 
  // part of the i-th sample if in the i-th row of the matrix there is a '1' at the column
  // that corresponds to this species (otherwise there is a '0').

  // The last argument is an output iterator of the type std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the value of the measure for the sample that is described in the i-th row of the input matrix.
  template <class OutputIterator>
  int matrix_query_weighted_basic(std::vector<std::string> &names, 
                                  std::vector< std::vector<bool> > &matrix, OutputIterator ot )
  {return this->_matrix_query_weighted(*p_tree, names, matrix, *this, false, ot);}


  // Same as the function above except that the returned values distance values have been "standardised":
  // from each value we have subtracted the mean value of the measure among all samples of leaves of the
  // same size, and we have divided the result by the deviation.
  template <class OutputIterator>
  int matrix_query_weighted_standardised(std::vector<std::string> &names, 
                                         std::vector< std::vector<bool> > &matrix, 
                                         OutputIterator ot, int repetitions=1000 )
  {
    if(p_tree->is_ultrametric() == true || this->probability_distribution() != Kernel::UNIFORM_FIXED_SIZE)
      return this->_matrix_query_weighted(*p_tree, names, matrix, *this, true, ot,repetitions);
    else 
      return -1;
  }


  // Input: A csv file which stores a matrix where each column corresponds to a species of the tree
  // and each row indicates a sample of these species for which we want to compute the
  // distance measure. A certain species is considered as part of the i-th sample if in the i-th row 
  // of the matrix there is a '1' at the column that corresponds to this species (otherwise there is a '0').

  // The second argument is an output iterator of the type std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the value of the measure for the sample that is described in the i-th row of the input matrix.
  template <class OutputIterator>
  int csv_matrix_query_weighted_basic(char *filename, OutputIterator ot )
  {return this->_csv_matrix_query_weighted(*p_tree, filename, *this, false, ot);}


  // Same as the function above except that the returned values distance values have been "standardised":
  // from each value we have subtracted the mean value of the measure among all samples of leaves of the
  // same size, and we have divided the result by the deviation.
  template <class OutputIterator>
  int csv_matrix_query_weighted_standardised(char *filename, OutputIterator ot, int repetitions=1000 )
  {
    if(p_tree->is_ultrametric() == true || this->probability_distribution() != Kernel::UNIFORM_FIXED_SIZE)
      return this->_csv_matrix_query_weighted(*p_tree, filename, *this, true, ot, repetitions);
    else 
      return -1;
  }

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


  // Computes all together the probability values f(x) = \binom{x}{r}/\binom{s}{r}

  void compute_all_hypergeometric_probabilities( int sample_size, int number_of_leaves);


  // Simulates f(x) = \binom{x}{r}/\binom{s}{r}, x \in [ _sample_size , _number_of_leaves ]

  Number_type hypergeom( int x )
  {
     if( x < _sample_size || x > _number_of_leaves )
       return Number_type(0.0);

     if( x == _number_of_leaves )
       return Number_type(1.0);

     return _hypergeom[x-_sample_size];
  }

  Number_type hypergeom_minus_one( int x )
  {
     if( _sample_size - 1 < x )
       return hypergeom(x)*Number_type(_sample_size)/Number_type(x-_sample_size + 1);

     if( _sample_size - 1 == x )
       return hypergeom(_sample_size);

     return Number_type(0.0);
  }

  Number_type hypergeom_minus_two( int x )
  {
     if( x < _sample_size-2)
       return Number_type(0.0);

     if( x == _sample_size-2)
       return hypergeom(_sample_size);

     if( x == _sample_size-1)
       return Number_type(_sample_size-1)*hypergeom(_sample_size);

     return hypergeom(x)*(Number_type(_sample_size)*
              Number_type(_sample_size-1))/(Number_type(x-_sample_size+1)*
              Number_type(x-_sample_size+2));
  }

  Number_type compute_expectation_uniform_distribution( int sample_size );

  Number_type compute_expectation( int sample_size );

  Number_type compute_variance_uniform_distribution
  ( int sample_size, Number_type expect = Number_type(-1.0) );
  
  Number_type compute_variance_uniform_distribution_slow
  ( int sample_size, Number_type expect = Number_type(-1.0) );

  Number_type compute_variance( int sample_size, Number_type expect = Number_type(-1.0) );

  Number_type compute_deviation( int sample_size, Number_type expect = Number_type(-1.0) );

  Number_type edge_contribution(int index, int sample_size_in_subtree, int sample_size)
  {
    if(!p_tree->is_ultrametric())
    {
      std::string exception_msg;
      exception_msg += " Request to compute edge contribution on a non-ultrametric tree.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    Node_type v = p_tree->node(index);

    if( sample_size_in_subtree < 0 || sample_size_in_subtree > v->all_subtree_leaves )
    {
      std::string exception_msg;
      exception_msg += 
           " Request to compute edge contribution with subtree sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(p_tree->is_root(index) || sample_size_in_subtree==0)
      return Number_type(0.0);

    if(sample_size_in_subtree>= sample_size)
      return Number_type(0.0);

    if(sample_size_in_subtree == 1)   
      return Number_type(2.0)*Number_type(v.distance);

    return Number_type(0.0);

  } // Number_type edge_contribution(int index, int sample_size_in_subtree)

  
  private:

  Tree_type                *p_tree; // Stores a pointer to a phylogenetic tree object.
  std::vector<Number_type> _hypergeom; // Stores f(x) = \binom{x}{r}/\binom{s}{r}, 
                                       // x \in [ _sample_size , _number_of_leaves ]

  std::vector<Number_type> _CPoisson_exps, _CPoisson_vars,
                           _Sequential_exps, _Sequential_devs;

  std::vector<Number_type> _max_subtree_path_costs; // The i-th slot of this vector stores
                                                    // the cost of the longest simple path between
                                                    // the node with index i and any marked leaf
                                                    // node in its subtree.

  int          _sample_size;
  int          _number_of_leaves;

}; // struct Mean_nearest_taxon_distance


} // namespace PhylogeneneticMeasures

#include "Mean_nearest_taxon_distance_impl.h"

#endif // MEAN_NEAREST_TAXON_DISTANCE_H
