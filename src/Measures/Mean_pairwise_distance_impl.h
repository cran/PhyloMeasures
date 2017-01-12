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

#ifndef MEAN_PAIRWISE_DISTANCE_IMPL_H
#define MEAN_PAIRWISE_DISTANCE_IMPL_H

#include<vector>

  template < class KernelType>
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  operator()( RangeIterator rbegin, RangeIterator rend,
              int min_index, int max_index )
  {
    if(rend-rbegin < 2)
      return Number_type(0.0);

    if(p_tree->number_of_nodes()<2)
      return Number_type(0.0);

    int intersection_index = p_tree->compute_intersection_node_index(min_index, max_index);

    // If the two paths coincide (intersection_index designates a leaf node)
    // then return zero distance.
    if( p_tree->node(intersection_index).children.size() == 0 )
      return Number_type(0.0);

    p_tree->node(intersection_index).mark = true;

    Number_type total_distance(0.0);
    int count_sample_nodes= int(rend-rbegin);

    // Mark all the nodes which fall on a path that connects
    // a leaf node in the sample and the node indicated by intersection_index.

    p_tree->mark_Steiner_tree_of_sample(rbegin, rend);
    p_tree->assign_marked_subtree_leaves(intersection_index);


    // Trace all edges in the marked subtree. Any edge that connects
    // a parent node v to child node u in this subtree appears as many
    // times in the final solution as k * (n-k) where k = u.marked_subtree_leaves
    // and n is the total number of nodes in the marked subtree.

    for( int i=1; i < p_tree->number_of_marked_nodes(); i++)
    {
      Node_type v = p_tree->node(p_tree->marked_node(i));

      total_distance += Number_type(v.distance)*
                        Number_type(v.marked_subtree_leaves)*
                        Number_type(count_sample_nodes - v.marked_subtree_leaves); 
    }

    // Unmark all marked nodes and their 'marked_subtree_leaves' fields.

    p_tree->unmark_Steiner_tree_of_sample(rbegin, rend);

    return total_distance*Number_type(2.0)/
           ( Number_type( count_sample_nodes)*Number_type(count_sample_nodes-1) );

  } // Mean_pairwise_distance<KernelType>::operator()(...)


  template< class KernelType >
  template < class OutputIterator >
  void PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  incremental_operator( std::vector<int> &sample,
                        std::vector<int> &sample_sizes, OutputIterator ot )
  {
    int n_l=p_tree->number_of_leaves(),
        n_n=p_tree->number_of_nodes();  

    for(int i=0; i<sample_sizes.size(); i++)
      if(sample_sizes[i] > n_l || sample_sizes[i] < 0 || sample_sizes[i] > sample.size())
      {
        std::string exception_msg;
        exception_msg += " One of the sample sizes provided to the incremental operator is out of range.\n";     
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }
      else if(i>0 && sample_sizes[i]<=sample_sizes[i-1])
      {
        std::string exception_msg;
        exception_msg += " The sample sizes provided to the incremental operator are not sorted.\n";     
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }

    if(sample_sizes.back() != sample.size())
    {
      std::string exception_msg;
      exception_msg += " The largest sample size should equal the size of the input sample.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    // Output values for all trivial sample sizes (i.e <2)
   
    int ss_index=0;

    while(ss_index < sample_sizes.size() && sample_sizes[ss_index] <2)
    {
      *ot++ = Number_type(0.0);
      ss_index++;
    }

    if(ss_index>=sample_sizes.size())
      return;

    if(sample_sizes.size()==0 || sample.size() == 0)
      return;

    // Deal with the first non-trivial sample size 

    this->initialize_marked_subtree_path_costs(*p_tree);

    int min_index=n_n+1, 
        max_index=-1, 
        prev_sample_size, 
        curr_sample_size;

    for(int i=0; i<sample_sizes[ss_index]; i++)
    {
      if(sample[i]<min_index)
        min_index = sample[i];

      if(sample[i]>max_index)
        max_index = sample[i];
    }

    int intersection_index = p_tree->compute_intersection_node_index(min_index, max_index);
    
    p_tree->node(intersection_index).mark = true;

    // Mark all the nodes which fall on a path that connects
    // a leaf node in the sample and the node indicated by intersection_index.
    // For each edge that we traverse, add the corresponding weight to the solution.

    std::vector<int>::iterator rit=sample.begin();

    for(int i=0; i<sample_sizes[ss_index]; i++)
      rit++;
      
    Number_type total_dist = 
    this->mark_tree_and_compute_subtree_path_costs(sample.begin(), rit, intersection_index);

    if( p_tree->node(intersection_index).children.size() > 0 )    
      *ot++ = Number_type(2.0)*total_dist/
              (Number_type(sample_sizes[ss_index])*Number_type(sample_sizes[ss_index]-1));
    else
      *ot++ = Number_type(0.0);

    // Continue with the rest sample sizes

    prev_sample_size = sample_sizes[ss_index];

    for(int i=ss_index+1; i<sample_sizes.size(); i++)
    {
      curr_sample_size = sample_sizes[i];

      for(int j=prev_sample_size; j<curr_sample_size; j++)
        total_dist += this->update_marked_subtree_path_costs(intersection_index,sample[j]); // intersection_index
                                                                                             // might get 
                                                                                             // updated.       

      *ot++ = Number_type(2.0)*total_dist/(Number_type(curr_sample_size)*Number_type(curr_sample_size-1));
  
      prev_sample_size = curr_sample_size;

    } // for(int i=ss_index+1; i<sample_sizes.size(); i++)

    // Unmark all marked nodes.

    this->clear_marked_subtree_path_costs(sample.begin(), sample.end());

    return;

  } // Mean_pairwise_distance<KernelType>::incremental_operator(...)

  // Consider the min cost Steiner tree that spans a given sample of leaf nodes (already marked). 
  // The following function computes for each node v of this subtree the sum of the costs of all paths
  // between v and the marked leaves in its subtree. This number is assigned to the slot of the vector
  // '_marked_subtree_path_costs' that has the same index as v. 
  // The field '.mark' of v is assumed to have already
  // been set to 'true' if the subtree of v contains at least one sample leaf.
  // Also, it is assumed thatthe fields 'marked_subtree_leaves' and 'marked_children'
  // have been already set to their correct values. The function returns a slightly 
  // different value than the '_marked_subtree_path_costs' of v; this is the sum of the costs of all 
  // paths that connect the father of v with the marked leaves in the subtree of v.

  template < class KernelType>
  typename KernelType::Number_type PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  assign_initial_marked_subtree_path_costs(int index, Number_type &mpd_dist)
  {  
    Node_type v = p_tree->node(index);

    // Check if v is a leaf node
    if( v.number_of_children() == 0 && index != p_tree->root_index() )
    {
      this->_marked_subtree_path_costs[index]=Number_type(0.0);
      return v.distance;
    }
    else if( v.number_of_children() == 0 )
    {
      this->_marked_subtree_path_costs[index]=Number_type(0.0);
      return Number_type(0.0);
    }
    else
    {
      // Interior node
      
      this->_marked_subtree_path_costs[index] = Number_type(0.0);
 
      for( int i = 0; i < v.number_of_marked_children(); i++ )
      {
        int curr_child = v.marked_children[i];
        Number_type child_cost = assign_initial_marked_subtree_path_costs(curr_child,mpd_dist);

        this->_marked_subtree_path_costs[index] += child_cost;
 
        mpd_dist += child_cost*Number_type(v.marked_subtree_leaves- p_tree->node(curr_child).marked_subtree_leaves);
      }   

    } // else of if( tree.node(index).number_of_children() == 0 )

    return this->_marked_subtree_path_costs[index] + (v.distance*Number_type(v.marked_subtree_leaves));
	  
  } // assign_initial_marked_subtree_path_costs(...)

  template < class KernelType>
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  mark_tree_and_compute_subtree_path_costs
  (RangeIterator rbegin, RangeIterator rend, int intersection_index)
  {
    p_tree->mark_Steiner_tree_of_sample(rbegin, rend);
    p_tree->assign_marked_subtree_leaves(intersection_index);

    Number_type mpd_dist(0.0);

    this->assign_initial_marked_subtree_path_costs(intersection_index, mpd_dist);

    return mpd_dist;
  }

  template < class KernelType>
  typename KernelType::Number_type PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  update_marked_subtree_path_costs(int &intersection_index, int new_node_index)
  {
    if(intersection_index == new_node_index)
      return Number_type(0.0);

    int old_intersection_index = intersection_index;
    int current_index, previous_index;
    Node_type previous_node;

    p_tree->update_marked_Steiner_tree(intersection_index, new_node_index);

    // If the new node is not a descendant of the intersection node,
    // update the cost data on the path between the intersection node and the
    // common ancestor of both nodes.

    if(old_intersection_index != intersection_index)
    {
      current_index = p_tree->node(old_intersection_index).parent;
      previous_index = old_intersection_index;
      previous_node = p_tree->node(previous_index);

      do
      {
        p_tree->insert_marked_node(current_index);
        p_tree->node(current_index).marked_subtree_leaves = previous_node.marked_subtree_leaves;

        this->_marked_subtree_path_costs[current_index]=
           this->_marked_subtree_path_costs[previous_index] +
           (previous_node.distance*previous_node.marked_subtree_leaves);

        previous_index = current_index;
        current_index = p_tree->node(current_index).parent;
        previous_node = p_tree->node(previous_index);

      }while(previous_index != intersection_index);

    } // if(old_intersection_index != intersection_index)

    Number_type path_dist(0.0),
                total_added_dist(0.0);

    current_index = p_tree->node(new_node_index).parent;
    previous_index = new_node_index; 
    previous_node = p_tree->node(new_node_index);

    p_tree->node(new_node_index).marked_subtree_leaves=1;
    p_tree->insert_marked_node(new_node_index);
    p_tree->node(new_node_index).mark = true;
    this->_marked_subtree_path_costs[new_node_index] = Number_type(0.0);

    // We have to insert to the marked_nodes vector the new nodes that appear 
    // in the tree and which are exclusive ancestors of the newly added leaf.
    // These fall inside the maximal chain of nodes that have degree exactly 
    // two and which connects the new leaf with rest of the tree.
    // Thus, the last of the nodes in this chain is connected with an edge
    // to a node with degree > 2. As soon as we find this condition, we stop
    // inserting nodes to the marked nodes vector.

    bool found_edge_to_existing_tree=false;

    do
    {
      // The next condition checks if the current index points to a 
      // node that has been freshly added to the marked subtree,
      // and it is different from the root of this subtree.

      if(found_edge_to_existing_tree==false && p_tree->node(current_index).number_of_marked_children()==1)
        p_tree->insert_marked_node(current_index);
      else if(found_edge_to_existing_tree==false)
        found_edge_to_existing_tree=true;
     
      path_dist += previous_node.distance;

      for(int i=0; i<p_tree->node(current_index).number_of_marked_children(); i++)
      {
        int curr_child = p_tree->node(current_index).marked_children[i];

        if(curr_child != previous_index)
          total_added_dist += this->_marked_subtree_path_costs[curr_child] +
                              ( p_tree->node(curr_child).distance*
                                Number_type(p_tree->node(curr_child).marked_subtree_leaves))+ 
                              (path_dist*p_tree->node(curr_child).marked_subtree_leaves);  
      }

      p_tree->node(current_index).marked_subtree_leaves++;
      this->_marked_subtree_path_costs[current_index]+= path_dist;

      previous_index = current_index;
      current_index = p_tree->node(current_index).parent;
      previous_node = p_tree->node(previous_index);

    }while(previous_index != intersection_index);

    return total_added_dist;

  } // update_marked_subtree_path_costs(...)

  template < class KernelType>
  template < class RangeIterator >
  void PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  clear_marked_subtree_path_costs(RangeIterator rbegin, RangeIterator rend)
  {
     for(int i=0; i<p_tree->number_of_marked_nodes(); i++)
       this->_marked_subtree_path_costs[p_tree->marked_node(i)] = Number_type(0.0);

     p_tree->unmark_Steiner_tree_of_sample(rbegin,rend);
  }

  // Computes for every leaf the average distance to all other leaves.
  // The output is a vector (interfaced by an output iterator) with
  // the average distance values. The first input argument is a vector 
  // with the leaf names indicating the order of the distances in the output.
  template < class KernelType>  
  template < class OutputIterator >
  void PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  compute_average_leaf_distances(std::vector<std::string> &leaf_names, OutputIterator  ot ) 
  {
    std::map<std::string, Number_type > names_to_distances;
      
    p_tree->assign_all_subtree_leaves();

    if(this->_edge_path_costs.size() == 0 )
      compute_all_edge_path_costs(*p_tree);

    for(Leaves_iterator lit = p_tree->leaves_begin(); lit != p_tree->leaves_end(); lit++)
      names_to_distances[lit->first]= this->_edge_path_costs[lit->second]; 

    for(int i=0; i<leaf_names.size(); i++)
      *ot++ = names_to_distances[leaf_names[i]]/Number_type(p_tree->number_of_leaves()-1);
  
  } // compute_average_leaf_distances(std::vector<std::string> &leaf_names, OutputIterator  ot )

  template < class KernelType>
  template < class OutputIterator >
  void PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  compute_average_leaf_distances_slow(std::vector<std::string> &leaf_names, OutputIterator  ot ) 
  {
    for(int i=0; i<leaf_names.size(); i++)
    {
      int index = p_tree->find_leaf(leaf_names[i])->second;
      Number_type dist(0.0);

      for(Leaves_iterator lit = p_tree->leaves_begin(); lit != p_tree->leaves_end(); lit++)      
        if(lit->second != index)
        {
          int index_2 = lit->second;
          std::vector<int> tmp_vec;
          tmp_vec.push_back(index);
          tmp_vec.push_back(index_2);
          dist += operator()(tmp_vec.begin(), tmp_vec.end(),std::min(index,index_2), std::max(index,index_2)); 
        }
            
      *ot++ = dist/Number_type(leaf_names.size()-1);
    }
  
  } // compute_average_leaf_distances_slow(std::vector<std::string> &leaf_names, OutputIterator  ot )

  template < class KernelType>
  typename KernelType::Number_type PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  compute_expectation_uniform_distribution( int sample_size )
  {
    if(sample_size < 0 || sample_size > p_tree->number_of_leaves())
    {
      std::string exception_msg;
      exception_msg += " Request to compute expectation with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(sample_size <= 1)
      return Number_type(0.0);

    if(_expectation_unif != Number_type(-1.0) )
      return _expectation_unif;

   int s(p_tree->number_of_leaves());

   _expectation_unif = Number_type(2.0)*(this->total_path_costs(*p_tree))/(Number_type(s)*Number_type(s-1));

    return _expectation_unif;
	
  } // compute_expectation_uniform_distribution( int sample_size )

  template < class KernelType>
  typename KernelType::Number_type PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  compute_variance_uniform_distribution( int sample_size )
  {
    if(sample_size < 0 || sample_size > p_tree->number_of_leaves())
    {
      std::string exception_msg;
      exception_msg += " Request to compute variance with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if( sample_size <= 1 || sample_size == p_tree->number_of_leaves() )
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(p_tree->root_index());

    int s(p_tree->number_of_leaves());

    Number_type c1 = Number_type(Number_type(4*(sample_size-2))*Number_type(sample_size-3))/
                                 (Number_type(s)*Number_type(sample_size)*
                                  Number_type(sample_size-1)*Number_type(s-1)*
                                  Number_type(s-2)*Number_type(s-3)),
                c2 = Number_type(4*(sample_size-2))/
                                  (Number_type(s)*Number_type(sample_size)*Number_type(sample_size-1)*Number_type(s-1)*
                                   Number_type(s-2)),
                c3 = Number_type(4)/ (Number_type(s)*Number_type(sample_size)*Number_type(sample_size-1)*
                                      Number_type(s-1));

    if( this->sum_all_edges_costs() == Number_type(-1.0) )
      this->compute_all_costs_values(*p_tree);

    return ( c1 * (this->total_path_costs(*p_tree))*(this->total_path_costs(*p_tree)) +
             (c2-c1)*(this->sum_all_leaf_costs()) +
             (c1-(Number_type(2.0)*c2)+c3)*(this->sum_all_edges_costs()) -
             compute_expectation(sample_size)*compute_expectation(sample_size) );
			 
  } // compute_variance_uniform_distribution( ... )

#endif //MEAN_PAIRWISE_DISTANCE_IMPL_H
