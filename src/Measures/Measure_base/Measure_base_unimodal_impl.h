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

#ifndef MEASURE_BASE_UNIMODAL_IMPL_H
#define MEASURE_BASE_UNIMODAL_IMPL_H

#include<string>
#include<vector>
#include<iostream>
#include<fstream>
#include<cstdlib>

  template < class KernelType >
  template < class TreeType, class OutputIterator >
  void PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _read_sample_sizes_from_file(char *filename, TreeType &tree, OutputIterator ot)
  {
    std::ifstream in(filename);
    std::vector<int> sample_sizes;

    // Reading first file with queries
    if( !( in.is_open() && in.good() ) )
    {
      std::string exception_msg(" There was a problem with opening the file with the species samples.\n");
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    std::string line;
    std::getline(in,line);
 
    int prev_index=-1, current_index =0;

    while( current_index< line.size()-1 )
    {
      do
      {
        current_index++;
      }
      while(line[current_index]!=',' && current_index < line.size() );

      if(current_index -prev_index < 2)
      {
        std::string exception_msg(" There is a mistake in the syntax of the sample sizes file.\n");
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      } 

      std::string substring = line.substr(prev_index+1,current_index -prev_index-1);

      for(int i=0; i<substring.size(); i++)
        if(isdigit(substring[i]) == false)
        {
          std::string exception_msg;
          exception_msg += " There is an error in the syntax of the sample sizes file.\n";     
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

      int size = atoi(substring.c_str());

      if(size < 0 || size > tree.number_of_leaves() )
      {
        std::string exception_msg(" One of the sample sizes in the file is out of range.\n");
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      } 

      sample_sizes.push_back(size);

      prev_index = current_index;

    } // while( current_index< line.size()-1 )

    for(int i=0; i<sample_sizes.size(); i++)
      *ot++ = sample_sizes[i];

  } // _read_sample_sizes_from_file(...)


  // Input:  A range of iterators that indicate a list of species names (in std::string format).
  // Output: The value of the current measure for this set of species.

  template < class KernelType >
  template < class TreeType, class RangeIterator, class Measure >    
  typename KernelType::Number_type 
  PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _list_query(TreeType &tree, RangeIterator rbegin, RangeIterator rend, Measure &msr)
  {
    typedef TreeType Tree_type;
    RangeIterator it;
    std::string str;
    typename Tree_type::Leaves_iterator lv_it;
    std::vector<int> leaf_indices;
    int min_index=tree.number_of_nodes(), max_index=-1;

    // Find the indices of all the leaf nodes that correspond to the query species.
    // Find also the minimum and maximum of those indices (here represented as min_index & max_index).
    for( it = rbegin; it != rend; it++ )
    {
      str = *it;
      lv_it = tree.find_leaf(str);
	  
      if( lv_it != tree.leaves_end() )
      {
         leaf_indices.push_back(lv_it->second);

         if( leaf_indices.back() < min_index )
           min_index = leaf_indices.back();

         if( leaf_indices.back() > max_index )
           max_index = leaf_indices.back();
      }
    }

    return msr(leaf_indices.begin(), leaf_indices.end(), min_index, max_index);
	
  } // list_query(TreeType &tree, RangeIterator rbegin, RangeIterator rend)
  
 
  // Input: A tree and a txt file that stores a list of species names, which constitute a subset 
  // of the species (leaf nodes) in the tree and which appear in random order.
  // Output: The value of the current measure for this set of species.
 
  template < class KernelType >
  template <class TreeType, class Measure>  
  typename KernelType::Number_type 
  PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _list_query(TreeType &tree, char* filename, Measure &msr)
  {
    std::vector<std::string> vec;

    std::ifstream in(filename);

    if( !( in.is_open() && in.good() ) )
    {
      std::string exception_msg(" There was a problem with opening the file with the species names.\n");
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    // Read tree string

    std::string name;
    char a;

    in >> a;

    while(in.good())
    {
      if( a == ',')
      {
        vec.push_back(name);
        name.clear();
      }
      else if( a != '\n' && a!='\0' )
        name.push_back(a);

      in >> a;
    }

    vec.push_back(name);

    in.close();

    return _list_query(tree,vec.begin(),vec.end(),msr);
  
  } // list_query(TreeType &tree, char* filename)

  // Input: a file name that stores a presence/absence csv matrix.
  // Output: a vector of std::strings storing the column names of the matrix,
  // and a 2D boolean matrix storing the actual matrix entries. 

  template < class KernelType >  
  void PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  read_csv_matrix(const char *filename, std::vector<std::string> &names, 
                  std::vector< std::vector<bool> > &matrix)
  {
    std::ifstream in(filename);

    if( !( in.is_open() && in.good() ) )
    {
      std::string exception_msg(" There was a problem with opening the file with the matrix.\n");
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }    

    names.clear();

    for(int i=0; i<matrix.size(); i++)
      matrix[i].clear();
  
    matrix.clear();

    // Read the first row, the one that contains the species names

    std::string line;
    char a;

    std::getline(in,line);

    int c=0;

    while(c < line.size() )
    {
      std::string str;
      a = line[c];

      while(a != ',' && a != ' ' && a != '\r' && c<line.size() )
      {
        str.push_back(a);
        c++;
        a = line[c];
      }

      if( str.size() > 0 )
        names.push_back(str);

      c++;

    } // while(c < line.size() )

    std::vector<bool> initialized_row;

    initialized_row.assign(names.size(),false);

    while( in.good() )
    {
      line.clear();

      int count=0, start=0;

      std::getline(in,line);

      if(!in.good())
        break;

      matrix.push_back(initialized_row);

      // Exclude the first word of the line if first character is not zero or one.

      if( line[start] != '0' && line[start] != '1' )
        while( start < line.size() && line[start] != ',' )
          start++;

      for( int i=start; i<line.size(); i++)
      {
        a = line[i];

        if( a == '1' )
        {
          if(count >= names.size())
          {
            std::string exception_msg;
            exception_msg += " The matrix file has wrong syntax.\n";
            Exception_type excp;
            excp.get_error_message(exception_msg);
            Exception_functor excf;
            excf(excp);
          }

          matrix.back()[count] = true;

          count++;
        }
        else if ( a == '0' )
        {
          if(count >= names.size())
          {
            std::string exception_msg;
            exception_msg += " The matrix file has wrong syntax.\n";
            Exception_type excp;
            excp.get_error_message(exception_msg);
            Exception_functor excf;
            excf(excp);
          }

          count++;
        }
        else if( a != ',' && a !=' ' && a !='\r' && a !='\n')
        {
          std::string exception_msg;
          exception_msg += " The matrix file has wrong syntax.\n";
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

      } // for( int i=start; i<line.size(); i++)

    } // while( in.good() )

    in.close();

  } // read_csv_matrix(...)


  // Input: A phylogenetic tree, a vector of std::strings storing species names, a boolean 2D matrix, and a 
  // distance measure. The j-th column of the matrix corresponds to a species in the tree 
  // (its name stored in the j-th entry of the names vector), and each row indicates a sample 
  // of these species for we want to compute the distance measure. A certain species is considered 
  // as part of the i-th sample if in the i-th row of the matrix there is a '1' at the column 
  // that corresponds to this species (otherwise a '0').

  // The last two arguments are boolean that indicates if it is required to compute the standardized
  // value of the measure, and an output iterator of the type:
  // std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the value of the measure for the sample that is described in the i-th row of the input matrix.
  
  template < class KernelType >  
  template < class TreeType, class Measure, class OutputIterator>
  int PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _matrix_query( TreeType &tree, std::vector<std::string> &names, 
                 std::vector< std::vector<bool> > &matrix, 
                 Measure &msr, bool standardised, OutputIterator ot )
  {
    typedef TreeType Tree_type;
    std::vector<Number_type> mean_values, deviation_values;
    std::vector<int> column_to_node_vec;

    mean_values.assign(tree.number_of_leaves()+1, Number_type(-1.0));
    deviation_values.assign(tree.number_of_leaves()+1, Number_type(-1.0));

    if(names.size() < tree.number_of_leaves())
    {
      std::string warning(" Warning: the input matrix has fewer columns than the number of species in the tree.");
      Exception_functor().issue_warning(warning);
    }

    std::vector<bool> checked_names;
    checked_names.assign(tree.number_of_nodes(),false);
 
    for( int i=0; i<names.size(); i++ )
    {
      typename Tree_type::Leaves_iterator lv_it = tree.find_leaf(names[i]);

      if( lv_it == tree.leaves_end() )
      {
        std::string exception_msg;
        exception_msg += " One of the species names in input the matrix was not found in the tree (";
        exception_msg += names[i];
        exception_msg += ") \n";  
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }
      else 
      {
        if(checked_names[(*lv_it).second] == true)
        {
          std::string exception_msg;
          exception_msg += " Two or more columns of the input matrix share the same species name (";
          exception_msg += (*lv_it).first;
          exception_msg += ") \n";  
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        } 
        else
          checked_names[(*lv_it).second] = true;
 
      } // else of if( lv_it == tree.leaves_end() )

      column_to_node_vec.push_back((*lv_it).second);
    }

    // Read the rest of the matrix, executing a query per line.

    for(int i=0; i<matrix.size(); i++)
    {
      std::vector<int> query_nodes; 
      int min = tree.number_of_nodes(), max = -1;

      for(int j=0; j<matrix[i].size(); j++)
        if(matrix[i][j] == true) 
        {
          query_nodes.push_back( column_to_node_vec[j] );

          if( query_nodes.back() < min )
            min = query_nodes.back();

          if( query_nodes.back() > max )
            max = query_nodes.back();
        }

      if( query_nodes.size() < 1 )
        *ot++ = 0.0;
      else
      {
        Number_type single_sample_result = msr(query_nodes.begin(), query_nodes.end(), min, max );

        if(standardised == false)  
          *ot++ = single_sample_result;
        else
        {
          if(mean_values[query_nodes.size()] == Number_type(-1.0))
          {
            mean_values[query_nodes.size()] = msr.compute_expectation(query_nodes.size());
            deviation_values[query_nodes.size()] = msr.compute_deviation(query_nodes.size());
          }  

          if(deviation_values[query_nodes.size()]==Number_type(0.0))
            *ot++ = single_sample_result - mean_values[query_nodes.size()];
          else
            *ot++ = (single_sample_result - mean_values[query_nodes.size()])/deviation_values[query_nodes.size()];

        } // else of if(standardised == false)
          
      } // else of if( query_nodes.size() < 1 )

    } // for(int i=0; i<matrix.size(); i++)

    return matrix.size();

  } // _matrix_query( ... )


  // Input: A csv file which stores a matrix where each column corresponds to a species of the tree
  // and each row indicates a sample of these species for we want to compute the
  // distance measure. A certain species is considered as part of the i-th sample if in the i-th row 
  // of the matrix there is a '1' at the column that corresponds to this species (otherwise a '0').

  // The second argument is an output iterator of the type: 
  //  std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the value of the measure for the sample that is described in the i-th row of the input matrix.
  
  template < class KernelType >  
  template < class TreeType, class Measure, class OutputIterator>
  int PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _csv_matrix_query( TreeType &tree, char *filename, Measure &msr, 
                     bool standardised, OutputIterator ot )
  {
    std::vector<std::string> names;
    std::vector<std::vector<bool> > matrix;

    this->read_csv_matrix(filename,names,matrix);
    return _matrix_query(tree, names, matrix, msr, standardised, ot);

  } // csv_matrix_query( ... )


  // Input: A phylogenetic tree, a vector of std::strings storing species names, a boolean 2D matrix, and a 
  // distance measure. The j-th column of the matrix corresponds to a species in the tree 
  // (its name stored in the j-th entry of the names vector), and each row indicates a sample 
  // of these species for we want to compute the distance measure. A certain species is considered 
  // as part of the i-th sample if in the i-th row of the matrix there is a '1' at the column 
  // that corresponds to this species (otherwise a '0').

  // The last two arguments are boolean that indicates if it is required to compute the standardized
  // value of the measure, and an output iterator of the type:
  // std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the value of the measure for the sample that is described in the i-th row of the input matrix.
  // If the 'standardised' argument is set to true then the returned values are  standardised values
  // of the measure for each row sample according to the Poisson binomial fixed-size distribution.


  // Precondition I: the leaves of the input tree store the probability values which are needed for
  // this distribution.

  // Precondition II: the probability_distribution of the input measure object is set to
  // Kernel::POISSON_BINOMIAL_FIXED_SIZE. 


  template < class KernelType >  
  template < class TreeType, class Measure, class OutputIterator>
  int PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _matrix_query_Poisson_binomial_fixed_size( TreeType &tree, std::vector<std::string> &names, 
                                             std::vector< std::vector<bool> > &matrix, 
                                             Measure &msr, bool standardised, OutputIterator ot )
  {
    if(!tree.stores_probability_values())
    {
      std::string exception_msg;
      exception_msg += " The leaves of the input tree do not store any probability values."; 
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    } 

    if(msr.probability_distribution() != Kernel::POISSON_BINOMIAL_FIXED_SIZE)
    {
      std::string exception_msg;
      exception_msg += " The distribution of the input measure object should be set to"; 
      exception_msg += " Kernel::POISSON_BINOMIAL_FIXED_SIZE .";
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    } 

    if(standardised == true)
    {
      // Scan the matrix and find which is the largest 
      // sample size (number of true entries) stored in its rows.
   
      int max_sum=0;

      for(int i=0; i<matrix.size(); i++)
      {
        int sum=0;

        for(int j=0; j<matrix[i].size(); j++)
          if(matrix[i][j]==true)
            sum++;

        if(i==0 || sum>max_sum)
          max_sum = sum;
 
      } // for(int i=0; i<matrix.size(); i++)

      // Compute the mean and variance of the input measure for the largest sample size
      // observed in the matrix. This will automatically compute also the mean and variances
      // for all smaller sample sizes.

      msr.compute_expectation(max_sum);

    } // if(standardised == true)

    return _matrix_query(tree, names, matrix, msr, standardised, ot);
    
  } // _matrix_query_Poisson_binomial_fixed_size


  template < class KernelType >  
  template < class TreeType, class Measure, class OutputIterator>
  int PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _csv_matrix_query_Poisson_binomial_fixed_size( TreeType &tree, char *filename, Measure &msr, 
                                                 bool standardised, OutputIterator ot )
  {
    std::vector<std::string> names;
    std::vector<std::vector<bool> > matrix;

    this->read_csv_matrix(filename,names,matrix);

    return _matrix_query_Poisson_binomial_fixed_size(tree, names, matrix, msr, standardised, ot);

  } // csv_matrix_query_Poisson_binomial_fixed_size( ... )



  // Input: A phylogenetic tree, a vector of std::strings storing species names, a boolean 2D matrix, and a 
  // distance measure. The j-th column of the matrix corresponds to a species in the tree 
  // (its name stored in the j-th entry of the names vector), and each row indicates a sample 
  // of these species for we want to compute the distance measure. A certain species is considered 
  // as part of the i-th sample if in the i-th row of the matrix there is a '1' at the column 
  // that corresponds to this species (otherwise a '0').

  // The last two arguments are: a boolean that indicates if it is required to compute the standardized
  // value of the measure, and an output iterator of the type:
  // std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the value of the measure for the sample that is described in the i-th row of the input matrix.
  // If the 'standardised' argument is set to true then the returned values are  standardised values
  // of the measure for each row sample according to the sequential fixed-size distribution.


  // Precondition I: the leaves of the input tree store the probability values which are needed for
  // this distribution.

  // Precondition II: the probability_distribution of the input measure object is set to
  // Kernel::SEQUENTIAL_FIXED_SIZE. 


  template < class KernelType >  
  template < class TreeType, class Measure, class OutputIterator>
  int PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _matrix_query_sequential_fixed_size( TreeType &tree, std::vector<std::string> &names, 
                                       std::vector< std::vector<bool> > &matrix, 
                                       Measure &msr, bool standardised, OutputIterator ot, int repetitions )
  {
    if(!tree.stores_probability_values())
    {
      std::string exception_msg;
      exception_msg += " The leaves of the input tree do not store any probability values."; 
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    } 

    if(msr.probability_distribution() != Kernel::SEQUENTIAL_FIXED_SIZE)
    {
      std::string exception_msg;
      exception_msg += " The distribution of the input measure object should be set to"; 
      exception_msg += " Kernel::SEQUENTIAL_FIXED_SIZE .";
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    } 

    if(standardised==false)
      return _matrix_query(tree, names, matrix, msr, standardised, ot);


    std::vector<int> leaf_indices;
    std::vector<Number_type> probabilities;
    typename TreeType::Leaves_iterator it;
 
    for(it = tree.leaves_begin(); it != tree.leaves_end(); it++)
    {
      leaf_indices.push_back(it->second);
      probabilities.push_back(tree.node_probability(it->second));
    }  

    Incremental_Monte_Carlo_handler  MC_handler;
    Sequential_sampler sampler(leaf_indices,probabilities);
    std::vector<std::pair<Number_type,Number_type> > moments;
    std::vector<Number_type > original_vals;

    MC_handler.estimate_moments_with_Monte_Carlo(msr, matrix, sampler, repetitions, std::back_inserter(moments));

    _matrix_query(tree, names, matrix, msr, false, std::back_inserter(original_vals));

    for(int i=0; i<original_vals.size(); i++)
    {
      if(moments[i].second==Number_type(0.0))
        *ot++ = original_vals[i] - moments[i].first;
      else
        *ot++ = (original_vals[i] - moments[i].first)/moments[i].second;
    }

    return matrix.size();

  } // _matrix_query_sequential_fixed_size


  // Input: A distance measure, a sample-size n,  and a number of repetitions. 

  // The other two arguments are two output iterators of the type:
  // std::back_insert_iterator<std::vector< Number_type > >.
  // Output: Two vector of n numbers each (passed in the form of the output iterators), where the i-th element
  // of the first vector is the expectation of the measure for sample-size i, and the i-th element of the
  // second vector is the deviation of the same measure.

  // Precondition I: the leaves of the measure's tree store the probability values which are needed for
  // this distribution.

  // Precondition II: the probability_distribution of the input measure object is set to
  // Kernel::SEQUENTIAL_FIXED_SIZE. 

  // Attention: this function might be too slow for lagre values of sample_size.
  // It is provided here only for completeness and for debugging reasons.

  template < class KernelType >  
  template < class Measure, class OutputIterator>
  void PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _compute_moments_sequential_fixed_size( Measure &msr, int sample_size, 
                                          OutputIterator ot_a, OutputIterator ot_b, int repetitions )
  {
    if(!msr.tree().stores_probability_values())
    {
      std::string exception_msg;
      exception_msg += " The leaves of the input tree do not store any probability values."; 
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    } 

    if(msr.probability_distribution() != Kernel::SEQUENTIAL_FIXED_SIZE)
    {
      std::string exception_msg;
      exception_msg += " The distribution of the input measure object should be set to"; 
      exception_msg += " Kernel::SEQUENTIAL_FIXED_SIZE .";
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    } 

    std::vector<int> leaf_indices;
    std::vector<Number_type> probabilities;
   
 
    typename Measure::Tree_type::Leaves_iterator it;
 
    for(it = msr.tree().leaves_begin(); it != msr.tree().leaves_end(); it++)
    {
      leaf_indices.push_back(it->second);
      probabilities.push_back(msr.tree().node_probability(it->second));
    }  

    Incremental_Monte_Carlo_handler  MC_handler;
    Sequential_sampler sampler(leaf_indices,probabilities);
    std::vector<std::pair<Number_type,Number_type> > moments;

    std::vector<int> query_sizes;
  
    for(int k=0; k<=sample_size; k++)
      query_sizes.push_back(k); 

    MC_handler.estimate_moments_with_Monte_Carlo(msr, query_sizes, sampler, 
                                                 repetitions, std::back_inserter(moments));

    for(int i=0; i<moments.size(); i++)
    {
      *ot_a++ = moments[i].first;
      *ot_b++ = moments[i].second;
    }

  } // _compute_moments_sequential_fixed_size

  template < class KernelType >  
  template < class TreeType, class Measure, class OutputIterator>
  int PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _csv_matrix_query_sequential_fixed_size( TreeType &tree, char *filename, Measure &msr, 
                                           bool standardised, OutputIterator ot, int repetitions)
  {
    std::vector<std::string> names;
    std::vector<std::vector<bool> > matrix;

    this->read_csv_matrix(filename,names,matrix);

    return _matrix_query_sequential_fixed_size(tree, names, matrix, msr, standardised, ot, repetitions);

  } // csv_matrix_query_sequential_fixed_size( ... )



  // Input: A phylogenetic tree, a vector of std::strings storing species names, a boolean 2D matrix, and a 
  // distance measure. The j-th column of the matrix corresponds to a species in the tree 
  // (its name stored in the j-th entry of the names vector), and each row indicates a sample 
  // of these species for we want to compute the distance measure. A certain species is considered 
  // as part of the i-th sample if in the i-th row of the matrix there is a '1' at the column 
  // that corresponds to this species (otherwise a '0').

  // The last argument is an output iterator of the type:
  // std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the p-value for the sample that is described in the i-th row of the input matrix,
  // computed according to the uniform fixed-size sampling distribution.


  // Precondition: the probability_distribution of the input measure object is set to
  // Kernel::UNIFORM_FIXED_SIZE. 



  template < class KernelType >  
  template < class TreeType, class Measure, class OutputIterator>
  int PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _pvalues_query_uniform_fixed_size( TreeType &tree, std::vector<std::string> &names, 
                                     std::vector< std::vector<bool> > &matrix, 
                                     Measure &msr, OutputIterator ot, int repetitions )
  {
    if(msr.probability_distribution() != Kernel::UNIFORM_FIXED_SIZE)
    {
      std::string exception_msg;
      exception_msg += " The distribution of the input measure object should be set to"; 
      exception_msg += " Kernel::UNIFORM_FIXED_SIZE .";
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    } 

    std::vector<int> leaf_indices;
    typename TreeType::Leaves_iterator it;
 
    for(it = tree.leaves_begin(); it != tree.leaves_end(); it++)
      leaf_indices.push_back(it->second);

    Incremental_Monte_Carlo_handler  MC_handler;
    Uniform_sampler sampler(leaf_indices);

    MC_handler.estimate_pvalues_with_Monte_Carlo(msr, names, matrix, sampler, repetitions, ot);
  
    return matrix.size();

  } // _pvalues_query_uniform_fixed_size


  template < class KernelType >  
  template < class TreeType, class Measure, class OutputIterator>
  int PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _csv_pvalues_query_uniform_fixed_size( TreeType &tree, const char *filename, Measure &msr, 
                                         OutputIterator ot, int repetitions)
  {
    std::vector<std::string> names;
    std::vector<std::vector<bool> > matrix;

    this->read_csv_matrix(filename,names,matrix);

    return _pvalues_query_uniform_fixed_size(tree, names, matrix, msr, ot, repetitions);

  } // _csv_pvalues_query_uniform_fixed_size( ... )


  // Input: A phylogenetic tree, a vector of std::strings storing species names, a boolean 2D matrix, and a 
  // distance measure. The j-th column of the matrix corresponds to a species in the tree 
  // (its name stored in the j-th entry of the names vector), and each row indicates a sample 
  // of these species for we want to compute the distance measure. A certain species is considered 
  // as part of the i-th sample if in the i-th row of the matrix there is a '1' at the column 
  // that corresponds to this species (otherwise a '0').

  // The last argument is an output iterator of the type:
  // std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the p-value for the sample that is described in the i-th row of the input matrix,
  // computed according to the R-style sequential weighted sampling distribution.


  // Precondition I: the leaves of the input tree store the probability values which are needed for
  // this distribution.

  // Precondition II: the probability_distribution of the input measure object is set to
  // Kernel::SEQUENTIAL_FIXED_SIZE. 


  template < class KernelType >  
  template < class TreeType, class Measure, class OutputIterator>
  int PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _pvalues_query_sequential_fixed_size( TreeType &tree, std::vector<std::string> &names, 
                                        std::vector< std::vector<bool> > &matrix, 
                                        Measure &msr, OutputIterator ot, int repetitions )
  {
    if(!tree.stores_probability_values())
    {
      std::string exception_msg;
      exception_msg += " The leaves of the input tree do not store any probability values."; 
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    } 

    if(msr.probability_distribution() != Kernel::SEQUENTIAL_FIXED_SIZE)
    {
      std::string exception_msg;
      exception_msg += " The distribution of the input measure object should be set to"; 
      exception_msg += " Kernel::SEQUENTIAL_FIXED_SIZE .";
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    } 

    std::vector<int> leaf_indices;
    std::vector<Number_type> probabilities;
    typename TreeType::Leaves_iterator it;
 
    for(it = tree.leaves_begin(); it != tree.leaves_end(); it++)
    {
      leaf_indices.push_back(it->second);
      probabilities.push_back(tree.node_probability(it->second));
    }  

    Incremental_Monte_Carlo_handler  MC_handler;
    Sequential_sampler sampler(leaf_indices,probabilities);

    MC_handler.estimate_pvalues_with_Monte_Carlo(msr, names, matrix, sampler, repetitions, ot);

    return matrix.size();

  } // _pvalues_query_sequential_fixed_size


  template < class KernelType >  
  template < class TreeType, class Measure, class OutputIterator>
  int PhylogeneticMeasures::Measure_base_unimodal<KernelType>::
  _csv_pvalues_query_sequential_fixed_size( TreeType &tree, const char *filename, Measure &msr, 
                                           OutputIterator ot, int repetitions)
  {
    std::vector<std::string> names;
    std::vector<std::vector<bool> > matrix;

    this->read_csv_matrix(filename,names,matrix);

    return _pvalues_query_sequential_fixed_size(tree, names, matrix, msr, ot, repetitions);

  } // _csv_pvalues_query_sequential_fixed_size( ... )

#endif //MEASURE_BASE_UNIMODAL_IMPL_H
