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

#ifndef INCREMENTAL_MONTE_CARLO_HANDLER_H
#define INCREMENTAL_MONTE_CARLO_HANDLER_H

#include<vector>
#include<set>
#include<algorithm>
#include<limits>
#include<chrono>
#include<random>
#include<thread>

namespace PhylogeneticMeasures
{
  
  template<class KernelType>
  class P_value_search_tree
  {

   public:

    typedef KernelType                               Kernel;
    typedef typename Kernel::Number_type             Number_type; 
    typedef typename Kernel::Numeric_traits          Numeric_traits;
    typedef typename Numeric_traits::Absolute_value  Absolute_value;

    struct Node_type 
    {
      Number_type max_val; 
      int lc, rc, parent;

      Node_type():max_val(-1.0), lc(-1),rc(-1),parent(-1){}

    }; // Node_type

   private:

    int _construct_tree()
    { 
      std::vector<int> nodes;
      int n=tree.size();
      int last_index=n-1;

      for(int i=0; i<n; i++)
        nodes.push_back(i);

      while(nodes.size() > 1)
      {
        std::vector<int> new_nodes;

        for(int i=0; i<nodes.size(); i++)
          if(i%2 == 1)
          {
            Node_type v;

            last_index++; 

            v.lc = nodes[i-1];
            v.rc = nodes[i];
            v.max_val = tree[nodes[i]].max_val;

            tree[nodes[i-1]].parent = last_index;
            tree[nodes[i]].parent = last_index;
            
            tree.push_back(v);

            new_nodes.push_back(last_index);
          }   

        if(nodes.size()%2 == 1)
          new_nodes.push_back(nodes.back());

       nodes = new_nodes;

      } // while(nodes.size() > 1)     

      return last_index; 

    } // _construct_tree(...)

    void _find_and_mark_recursive(Number_type value, int index)
    {
      Node_type v = tree[index];

      if(v.lc==-1)
      {
        if(v.max_val<=value || Absolute_value()(value-v.max_val) < 
         Number_type(0.01))//Number_type(10.0)*std::sqrt(std::numeric_limits<double>::epsilon()) )
         vdist[index+1].first++;
        else
         vdist[index].first++;       

        return;
      }

      if(tree[v.lc].max_val<=value || 
         Absolute_value()(value-tree[v.lc].max_val) < 
         Number_type(0.01))//Number_type(10.0)*std::sqrt(std::numeric_limits<double>::epsilon()) )
        _find_and_mark_recursive(value,v.rc);
      else
        _find_and_mark_recursive(value,v.lc);

    } // void _find_and_mark_recursive(Number_type value, int index)

   public:

    // The first element of each input pair is the value of the examined measure
    // for a given row of the input matrix. The second element is the index of
    // this row in the matrix. 
    P_value_search_tree(std::vector< std::pair<Number_type,int> > &vals, Number_type &ref_value)
    {

      std::pair<int,int> first_pr(0,-1);
      
      vdist.push_back(first_pr);

      for(int i=0; i<vals.size(); i++)
      {
        Node_type nd;
      
        nd.max_val = vals[i].first;
        tree.push_back(nd);
        std::pair<int,int> pr(0,vals[i].second);
        vdist.push_back(pr);
      }      
  
      _ref_value = ref_value;

      _construct_tree();

    } // P_value_search_tree(std::vector<Number_type> &vals)

    void find_and_mark(Number_type value)
    { _find_and_mark_recursive(value,tree.size()-1); }

   public:

    std::vector<Node_type> tree;
    std::vector< std::pair<int,int> > vdist; // First element is the partial rank (p-value), 
                                             // second element is the corresponding row in the
                                             // input matrix.
    Number_type _ref_value;

  }; // class P_value_search_tree

template<typename KernelType>
class Incremental_Monte_Carlo_handler
{
 public:

  typedef KernelType                                Kernel;
  typedef typename Kernel::Number_type              Number_type;
  typedef P_value_search_tree<Kernel>               Search_tree;
  typedef typename Kernel::Numeric_traits           Numeric_traits;
  typedef typename Numeric_traits::Absolute_value   Absolute_value;
  typedef typename Numeric_traits::Square_root      Square_root;

 protected:

  struct Is_smaller_pair
  {
    bool operator()(const std::pair<int,int>&a, const std::pair<int,int> &b)
    {
      if(a.first < b.first)
        return true;

      if(a.first > b.first)
        return false;

      if(a.second < b.second)
        return true;

      return false;
    }

  }; // struct Is_smaller_pair
  
  struct Has_smaller_value
  {
    bool operator()(const std::pair<Number_type,int>&a, 
                      const std::pair<Number_type,int> &b)
    {
      if(a.first < b.first)
        return true;

      if(a.first > b.first)
        return false;

      if(a.second < b.second)
        return true;

      return false;
    }

  }; // struct Has_smaller_value


  template< class Measure_Type, class Sampler_Type>
  struct Incremental_moments_functor
  {
   public:

    Incremental_moments_functor
    ( Measure_Type *msr, Sampler_Type *sampler, 
      std::vector<int> *sample_sizes, int repetitions, 
      std::vector<Number_type> *sums, 
      std::vector<Number_type> *square_sums):
      _msr(msr), _sampler(sampler), _sample_sizes(sample_sizes),  
      _sums(sums), _square_sums(square_sums), _repetitions(repetitions)
    {}
           
    void operator()(void) 
    {
      for(int i=0; i<_repetitions; i++)
      {
        std::vector<int> sample;
        std::vector<Number_type> res;

        (*_sampler)(_sample_sizes->back(), sample);

        _msr->incremental_operator(sample,*_sample_sizes, std::back_inserter(res));

        for(int j=0; j<res.size(); j++)
        {
          (*_sums)[j] += res[j];
          (*_square_sums)[j] += (res[j]*res[j]);
        }

      } // for(int i=0; i<_repetitions; i++)
  
    } // void operator()(void)

    Measure_Type *_msr; 
    Sampler_Type *_sampler;
    std::vector<int> *_sample_sizes; 
    std::vector<Number_type> *_sums; 
    std::vector<Number_type> *_square_sums; 
    int _repetitions;    

  }; // struct Incremental_moments_functor

  template< class Measure_Type, class Sampler_Type>
  struct Incremental_pvalues_functor
  {
   public:

    Incremental_pvalues_functor
    ( Measure_Type *msr, Sampler_Type *sampler, std::vector<int> *sample_sizes, 
      std::vector<Search_tree> *search_trees, int repetitions):
      _msr(msr), _sampler(sampler), _sample_sizes(sample_sizes), 
      _search_trees(search_trees), _repetitions(repetitions){}
           
    void operator()(void) 
    {
      for(int i=0; i<_repetitions; i++)
      {
        std::vector<int> sample;
        std::vector<Number_type> res;

        (*_sampler)(_sample_sizes->back(), sample);

        _msr->incremental_operator(sample,*_sample_sizes, std::back_inserter(res));

        for(int j=0; j<res.size(); j++)
          (*_search_trees)[j].find_and_mark(res[j]);

      } // for(int i=0; i<_repetitions; i++)

      // Aggregate the marked ranks over each search tree

      for(int i=0; i<(*_search_trees).size(); i++)
      {
        int total_count=0;

        for(int j=(*_search_trees)[i].vdist.size()-1; j>=0; j--)
        {
          (*_search_trees)[i].vdist[j].first+=total_count;
          total_count = (*_search_trees)[i].vdist[j].first;
        }    

      } // for(int i=0; i<(*_search_trees).size(); i++)

    } // void operator()(void)

    Measure_Type *_msr; 
    Sampler_Type *_sampler;
    std::vector<int> *_sample_sizes; 
    std::vector<Search_tree> *_search_trees; 
    int _repetitions; 

  }; // struct Incremental_pvalues_functor

 public:

  Incremental_Monte_Carlo_handler(){}

  template<class OutputIterator>
  void extract_sample_sizes(std::vector<int> &query_sizes, OutputIterator ot);

  template<class OutputIterator>
  void extract_sample_sizes(std::vector<std::vector<bool> > &matrix, OutputIterator ot);

  template<class OutputIterator>
  void extract_sample_sizes(std::vector<std::vector<bool> > &matrix_a,
                             std::vector<std::vector<bool> > &matrix_b,
                             bool is_symmetric, OutputIterator ot);

  template<class OutputIterator>
  void extract_sample_sizes_specific_pairs(std::vector<std::vector<bool> > &matrix_a,
                                           std::vector<std::vector<bool> > &matrix_b,
                                           std::vector<std::pair<int,int> > row_pairs,
                                           OutputIterator ot);
 
  template<class Measure>
  void extract_sample_size_sets
  (std::vector<std::string> &names, std::vector<std::vector<bool> > &matrix, 
   Measure &msr, int number_of_leaves, std::vector<int> &sample_sizes, 
   std::vector<std::vector<std::pair<Number_type, int> > > &vals_and_indices)
  {
    sample_sizes.clear();
    vals_and_indices.clear();

    std::vector<std::vector<std::pair<Number_type, int> > > range;
   
    range.assign(number_of_leaves+1,std::vector<std::pair<Number_type, int> >());

    std::vector<Number_type> vals;
    std::vector<int> sample_sizes_tmp;

    msr.matrix_query_basic(names, matrix, std::back_inserter(vals));


    for(int i=0; i<matrix.size(); i++)
    {
      int sample_size=0;

      for(int j=0; j<matrix[i].size(); j++)
        if(matrix[i][j]==true)
          sample_size++;

      range[sample_size].push_back(std::make_pair(vals[i],i));

    } // for(int i=0; i<matrix.size(); i++)

    
    for(int i=0; i<range.size(); i++)
      if(range[i].size() != 0)
      {
        sample_sizes.push_back(i);

        std::sort(range[i].begin(), range[i].end(), Has_smaller_value());

        vals_and_indices.push_back(range[i]);
      }    

    return;

  } // extract_sample_size_sets(...)

  template<class Measure, class SamplerType, class OutputIterator>
  void estimate_moments_with_Monte_Carlo
  (Measure &msr, std::vector<std::vector<bool> > &matrix, 
   SamplerType &sampler, int repetitions, OutputIterator ot)
  {
    std::vector<int> query_sizes;

    query_sizes.assign(matrix.size(), 0);

    for(int i=0; i<matrix.size(); i++)
      for(int j=0; j<matrix[i].size(); j++)
        if(matrix[i][j]==true)
          query_sizes[i]++;
 
    this->estimate_moments_with_Monte_Carlo(msr, query_sizes, sampler, repetitions, ot);   

  } // estimate_moments_with_Monte_Carlo( ... , std::vector<std::vector<bool> > &matrix, ... )

  template<class Measure, class SamplerType, class OutputIterator>
  void estimate_moments_with_Monte_Carlo
  (Measure &msr, std::vector<int> &query_sizes, 
   SamplerType &sampler, int repetitions, OutputIterator ot) 
  {
    unsigned seed;

    typedef typename SamplerType::Sampler_data  Sampler_data; 

    std::vector<int> sample_sizes;

    this->extract_sample_sizes(query_sizes, std::back_inserter(sample_sizes));

    int number_of_threads = std::max(int(std::thread::hardware_concurrency()),1);
    std::vector<std::thread> threads;

    std::vector< std::vector<Number_type> > sums, square_sums;

    sums.assign(number_of_threads, std::vector<Number_type>());

    // Matrices that will store the results
 
    for(int i=0; i<number_of_threads; i++)
      sums[i].assign(sample_sizes.size(), Number_type(0.0));

    square_sums.assign(number_of_threads, std::vector<Number_type>());

    for(int i=0; i<number_of_threads; i++)
      square_sums[i].assign(sample_sizes.size(), Number_type(0.0));

    // Create an integer generator which is meant only for producing seeds

    if(msr.seed()<0)
      seed = std::chrono::system_clock::now().time_since_epoch().count();    
    else
      seed = msr.seed();

    std::default_random_engine seed_generator(seed);
    std::uniform_int_distribution<unsigned int> seed_distribution(0,std::numeric_limits<unsigned int>::max());

    typedef Incremental_moments_functor<Measure,SamplerType>  Incremental_moments_functor;
    typedef typename Measure::Tree_type                      Tree_type;

    std::vector<Tree_type> tree_copies;
    std::vector<Measure> measure_copies;
    std::vector<SamplerType> sampler_copies;

    for(int i=0; i<number_of_threads; i++ )
    {
      Tree_type tree_tmp = msr.tree();
      tree_copies.push_back(tree_tmp);
    }

    for(int i=0; i<number_of_threads; i++ )
    {
      Measure msr_tmp(tree_copies[i]); 

      msr_tmp.set_auxiliary_data(msr.auxiliary_data());

      measure_copies.push_back(msr_tmp);

      unsigned int thr_seed = seed_distribution(seed_generator);

      SamplerType thr_sampler(sampler.data(),thr_seed);

      sampler_copies.push_back(thr_sampler);

    } // for(int i=0; i<number_of_threads; i++ )

    for(int i=0; i<number_of_threads; i++ )
    {
      int thr_repetitions = repetitions/number_of_threads;

      if(i<repetitions%number_of_threads)
        thr_repetitions++;

      Incremental_moments_functor imf(&(measure_copies[i]), &(sampler_copies[i]), &sample_sizes, 
                                      thr_repetitions, &(sums[i]), &(square_sums[i]));

      threads.push_back(std::thread(imf));

    }  // for(int i=0; i<number_of_threads; i++ )


    for(int i=0; i<threads.size(); i++)
      threads[i].join();

    threads.clear();

    std::vector<Number_type> means, variances;

    means.assign(sample_sizes.size(), Number_type(0.0));
    variances.assign(sample_sizes.size(), Number_type(0.0));

    for(int i=0; i<number_of_threads; i++)
      for(int j=0; j<sums[i].size(); j++)
      {
        means[j]= means[j] + sums[i][j]; 
        variances[j]= variances[j] + square_sums[i][j];   
      }

    for(int i=0; i<sample_sizes.size(); i++)
    {
      means[i] = means[i]/Number_type(repetitions); 
      variances[i] = (Number_type(repetitions)/Number_type(repetitions-1))*
                     ((variances[i]/Number_type(repetitions)) - (means[i]*means[i])); 

      if(variances[i] < Number_type(0.0))
        variances[i] = Number_type(0.0);
    }

    std::vector<int> size_to_index;

    size_to_index.assign(msr.tree().number_of_leaves()+1, -1);

    for(int i=0; i<sample_sizes.size(); i++)
      size_to_index[sample_sizes[i]] = i;

    for(int i=0; i<query_sizes.size(); i++)
      *ot++ = std::make_pair(means[size_to_index[query_sizes[i]]],
                             Square_root()(variances[size_to_index[query_sizes[i]]]));

    for(int i=0; i<tree_copies.size(); i++)
      tree_copies[i].clear();

  } // estimate_moments_with_Monte_Carlo(...)

  // Definition of p-value for observed value v compared to array of randomised values R[1 ... n]:
  // pval(v) = (1+ "number of values in R that are larger than or equal to v")/(n+1)
   
  template<class Measure, class SamplerType, class OutputIterator>
  void estimate_pvalues_with_Monte_Carlo
  (Measure &msr, std::vector<std::string> &names, 
   std::vector<std::vector<bool> > &matrix, 
   SamplerType &sampler, int repetitions, OutputIterator ot) 
  {
    typedef typename SamplerType::Sampler_data  Sampler_data; 

    std::vector<int> sample_sizes;
    std::vector<std::vector<std::pair<Number_type, int> > > vals_and_indices;
    unsigned seed;

    this->extract_sample_size_sets(names, matrix, msr, msr.tree().number_of_leaves(), 
                                   sample_sizes, vals_and_indices);

    int number_of_threads = std::max(int(std::thread::hardware_concurrency()),1);
    std::vector<std::thread> threads;

    Number_type ref_value = msr.reference_value();

    // Construct vectors of search trees. Each vector is to 
    // be assigned to a different processor, each search tree 
    // in a vector corresponding to a different sample size s
    // (it stores the values of the measure for the rows in
    //  the input matrix that have sample size s).

    std::vector< std::vector<Search_tree> > search_tree_copies;

    std::vector<Search_tree> vec_stree;

    for(int i=0; i<sample_sizes.size(); i++)
    {
      Search_tree stree(vals_and_indices[i], ref_value);
      vec_stree.push_back(stree);
    }

    for(int k=0; k<number_of_threads; k++)
      search_tree_copies.push_back(vec_stree);

    // Create an integer generator which is meant only for producing seeds

    if(msr.seed()<0)
      seed = std::chrono::system_clock::now().time_since_epoch().count(); 
    else
      seed = msr.seed();

    std::default_random_engine seed_generator(seed);
    std::uniform_int_distribution<unsigned int> seed_distribution(0,std::numeric_limits<unsigned int>::max());

    typedef Incremental_pvalues_functor<Measure,SamplerType>  Incremental_pvalues_functor;
    typedef typename Measure::Tree_type  Tree_type;

    std::vector<Tree_type> tree_copies;   
    std::vector<Measure> measure_copies;
    std::vector<SamplerType> sampler_copies;

    for(int i=0; i<number_of_threads; i++ )
    {
      Tree_type tree_tmp(msr.tree());
      tree_copies.push_back(tree_tmp);
    }

    for(int i=0; i<number_of_threads; i++ )
    {
      Measure msr_tmp(tree_copies[i]); 

      msr_tmp.set_auxiliary_data(msr.auxiliary_data());

      measure_copies.push_back(msr_tmp);

      unsigned int thr_seed = seed_distribution(seed_generator);

      SamplerType thr_sampler(sampler.data(),thr_seed);
      sampler_copies.push_back(thr_sampler);

    } // for(int i=0; i<number_of_threads; i++ )


    for(int i=0; i<number_of_threads; i++ )
    {
      int thr_repetitions = repetitions/number_of_threads;

      if(i<repetitions%number_of_threads)
        thr_repetitions++;

      Incremental_pvalues_functor ipf(&(measure_copies[i]), &(sampler_copies[i]), &sample_sizes,
                                      &(search_tree_copies[i]), thr_repetitions);

      threads.push_back(std::thread(ipf));

    }  // for(int i=0; i<number_of_threads; i++ )

    for(int i=0; i<threads.size(); i++)
      threads[i].join();

    threads.clear();

    std::vector<int> pvalue_ranks;

    pvalue_ranks.assign(matrix.size(), 0);

    for(int i=0; i<search_tree_copies.size(); i++)
      for(int j=0; j<search_tree_copies[i].size(); j++)
        for(int k=1; k<search_tree_copies[i][j].vdist.size(); k++)
          pvalue_ranks[search_tree_copies[i][j].vdist[k].second]+=search_tree_copies[i][j].vdist[k].first;   

    for(int i=0; i<pvalue_ranks.size(); i++)
      *ot++ = Number_type(pvalue_ranks[i]+1)/Number_type(repetitions+1);

    for(int i=0; i<tree_copies.size(); i++)
      tree_copies[i].clear();

    search_tree_copies.clear();

    sampler_copies.clear();
   
  } // estimate_pvalues_with_Monte_Carlo(...)

}; // class Incremental_Monte_Carlo_handler

} // namespace PhylogeneticMeasures

#include "Incremental_Monte_Carlo_handler_impl.h"

#endif // INCREMENTAL_MONTE_CARLO_HANDLER_H
