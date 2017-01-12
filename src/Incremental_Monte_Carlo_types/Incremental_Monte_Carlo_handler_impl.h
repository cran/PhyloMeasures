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

#ifndef INCREMENTAL_MONTE_CARLO_HANDLER_IMPL_H
#define INCREMENTAL_MONTE_CARLO_HANDLER_IMPL_H

  template< class KernelType >
  template<class OutputIterator>
  void PhylogeneticMeasures::Incremental_Monte_Carlo_handler<KernelType>::
  extract_sample_sizes(std::vector<std::vector<bool> > &matrix, OutputIterator ot)
  {
    std::set<int> sums;

    for(int i=0; i<matrix.size(); i++)
    {
      int s=0;

      for(int j=0; j<matrix[i].size(); j++)
        if(matrix[i][j] == true )
          s++;

      sums.insert(s);
    }

    typename std::set<int>::iterator it;

    for(it=sums.begin(); it!=sums.end(); it++)
      *ot++ = *it;

  } // extract_sample_sizes(...) // single matrix

  template< class KernelType >
  template<class OutputIterator>
  void PhylogeneticMeasures::Incremental_Monte_Carlo_handler<KernelType>::
  extract_sample_sizes(std::vector<int> &query_sizes, OutputIterator ot)
  {
    std::set<int> sizes;

    for(int i=0; i<query_sizes.size(); i++)
      sizes.insert(query_sizes[i]);

    typename std::set<int>::iterator it;

    for(it=sizes.begin(); it!=sizes.end(); it++)
      *ot++ = *it;

  } // extract_sample_sizes(...) // single matrix

  template< class KernelType >
  template<class OutputIterator>
  void PhylogeneticMeasures::Incremental_Monte_Carlo_handler<KernelType>::
  extract_sample_sizes(std::vector<std::vector<bool> > &matrix_a,
                       std::vector<std::vector<bool> > &matrix_b,
                       bool is_symmetric, OutputIterator ot)
  {
    std::set<int> sums_a, sums_b;

    for(int i=0; i<matrix_a.size(); i++)
    {
      int s=0;

      for(int j=0; j<matrix_a[i].size(); j++)
        if(matrix_a[i][j] == true )
          s++;

      sums_a.insert(s);
    }

    for(int i=0; i<matrix_b.size(); i++)
    {
      int s=0;

      for(int j=0; j<matrix_b[i].size(); j++)
        if(matrix_b[i][j] == true )
          s++;

      sums_b.insert(s);
    }

    typename std::set<int>::iterator it_a, it_b;

    if(is_symmetric)
    {
      std::set< std::pair<int,int>, Is_smaller_pair > size_pairs;

      for(it_a=sums_a.begin(); it_a != sums_a.end(); it_a++)
        for(it_b=sums_b.begin(); it_b != sums_b.end(); it_b++)
        {
          std::pair<int,int> pr(*it_a, *it_b);

          if(pr.first > pr.second)
          {
            int tmp = pr.second;
            pr.second = pr.first;
            pr.first = tmp;
          }  
   
          size_pairs.insert(pr);

        } // for(it_b=sums_b.begin(); it_b != sums_b.end(); it_b++)


      typename std::set<std::pair<int,int>,Is_smaller_pair>::iterator it_pr;

      for(it_pr=size_pairs.begin(); it_pr!=size_pairs.end(); it_pr++)
        *ot++ = *it_pr;

    } // if(is_symmetric)
    else
      for(it_a=sums_a.begin(); it_a != sums_a.end(); it_a++)
        for(it_b=sums_b.begin(); it_b != sums_b.end(); it_b++)
          *ot++ = std::make_pair(*it_a, *it_b);

  } // extract_sample_sizes(...) // two matrices


  template< class KernelType >
  template<class OutputIterator>
  void PhylogeneticMeasures::Incremental_Monte_Carlo_handler<KernelType>::
  extract_sample_sizes_specific_pairs(std::vector<std::vector<bool> > &matrix_a,
                                      std::vector<std::vector<bool> > &matrix_b,
                                      std::vector<std::pair<int,int> > row_pairs,
                                      OutputIterator ot)
  {
    std::vector<int> sums_a, sums_b;
    std::set< std::pair<int,int>, Is_smaller_pair > size_pairs;

    for(int i=0; i<matrix_a.size(); i++)
    {
      int s=0;

      for(int j=0; j<matrix_a[i].size(); j++)
        if(matrix_a[i][j] == true )
          s++;

      sums_a.push_back(s);
    }

    for(int i=0; i<matrix_b.size(); i++)
    {
      int s=0;

      for(int j=0; j<matrix_b[i].size(); j++)
        if(matrix_b[i][j] == true )
          s++;

      sums_b.push_back(s);
    }
 
    for(int i=0; i<row_pairs.size(); i++)
    {
       std::pair<int,int> pr(sums_a[row_pairs[i].first],sums_b[row_pairs[i].second]);
       size_pairs.insert(pr);
    }  
 
    typename std::set<std::pair<int,int>,Is_smaller_pair>::iterator it_pr;

    for(it_pr=size_pairs.begin(); it_pr!=size_pairs.end(); it_pr++)
      *ot++ = *it_pr;

  } // extract_sample_sizes_specific_pairs(...)

#endif // INCREMENTAL_MONTE_CARLO_HANDLER_IMPL_H
