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

#ifndef MEAN_PAIRWISE_DISTANCE_BASE_H
#define MEAN_PAIRWISE_DISTANCE_BASE_H

#include<vector>

namespace PhylogeneticMeasures {

template <class KernelType, class TreeType>
class Mean_pairwise_distance_base
{
  typedef KernelType                                     Kernel;
  typedef TreeType                                       Tree_type;
  typedef Mean_pairwise_distance_base<Kernel,Tree_type>  Self;
  typedef typename Kernel::Number_type                  Number_type; 
  typedef typename Tree_type::Node_type                 Node_type;

 public:
  
   Mean_pairwise_distance_base():_total_path_costs(-1.0),_sum_all_edges_costs(-1.0),_sum_all_leaf_costs(-1.0)
   {}

   Self& operator=(const Self& d)
   {
     _edge_path_costs.clear();
     _marked_subtree_path_costs.clear();

     for(int i=0; i<d._edge_path_costs.size(); i++)
       _edge_path_costs.push_back(d.edge_path_cost(i));

     for(int i=0; i<d._marked_subtree_path_costs.size(); i++)
       _marked_subtree_path_costs.push_back(d.edge_path_cost(i));

     _total_path_costs=d._total_path_costs; 
     _sum_all_edges_costs=d.sum_all_edges_costs(); 
     _sum_all_leaf_costs=d.sum_all_leaf_costs();

     return *this;

   } // operator=(const Self& d)

  void initialize_marked_subtree_path_costs(Tree_type &tree)  
  {
    if(_marked_subtree_path_costs.size() != tree.number_of_nodes())
    {  
      _marked_subtree_path_costs.clear();
      _marked_subtree_path_costs.assign(tree.number_of_nodes(), Number_type(0.0));
    }
  }

  Number_type total_path_costs(Tree_type &tree);

  // Helper function that computes all the sums-of-path-costs
  // values that are stored for future use.
  void compute_all_costs_values(Tree_type& tree);

  void compute_all_edge_path_costs(Tree_type& tree);

  Number_type edge_path_cost( int i ) const
  { return _edge_path_costs[i]; }

  Number_type sum_all_edges_costs( void ) const 
  { return _sum_all_edges_costs; }

  Number_type sum_all_leaf_costs( void ) const
  { return _sum_all_leaf_costs; }  

 protected:

  Number_type
  _compute_single_edge_path_costs( Tree_type& tree, int node_index, Number_type sum_anc1, Number_type sum_anc2,
                                   Number_type& all_weights );

 protected:

   std::vector<Number_type>  _edge_path_costs; // The i-th element of this vector stores
                                               // the sum of the costs of all paths that
                                               // contain the edge/node with index i in the
                                               // tree on which this class was applied.

   std::vector<Number_type>  _marked_subtree_path_costs; // The i-th element of this vector stores
                                                         // the sum of the costs of all paths between
                                                         // the node with index i in the associated tree
                                                         // and the marked leaves in its subtree.


   // The following are constants that are stored so as to achieve
   // constant time execution for the computation of the standard deviation

   Number_type _total_path_costs, _sum_all_edges_costs, _sum_all_leaf_costs;

}; // class Mean_pairwise_distance_base

} // namespace PhylogeneticMeasures

#include "Mean_pairwise_distance_base_impl.h"

#endif // MEAN_PAIRWISE_DISTANCE_BASE_H
