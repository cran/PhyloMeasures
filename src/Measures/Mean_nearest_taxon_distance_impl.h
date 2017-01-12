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

#ifndef MEAN_NEAREST_TAXON_DISTANCE_IMPL_H
#define MEAN_NEAREST_TAXON_DISTANCE_IMPL_H

  // The following function computes the two smallest costs of paths
  // that connect the node with index `current_index' to a leaf in its subtree.
  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  _compute_subtree_min_values( Tree_type &tree, int current_index )
  {
    Node_type current_node = tree.node(current_index);

    for( int i=0; i<current_node.marked_children.size(); i++ )
    {
      Number_type cmin = _compute_subtree_min_values(tree, current_node.marked_children[i]);

      if( tree.node(current_index).first_min == Number_type(-1.0) || cmin < tree.node(current_index).first_min)
      {
        tree.node(current_index).second_min = tree.node(current_index).first_min;
        tree.node(current_index).first_min = cmin;
      }
      else if(tree.node(current_index).second_min == Number_type(-1.0) || 
              cmin < tree.node(current_index).second_min )
        tree.node(current_index).second_min = cmin;

    } // for( int i=0; i<current_node.marked_children.size(); i++ )

    if( tree.node(current_index).marked_children.size() == 0 )
    {
      tree.node(current_index).first_min = Number_type(0.0);
      tree.node(current_index).second_min = Number_type(0.0);
    }

    return tree.node(current_index).first_min + Number_type(tree.node(current_index).distance);

  } // _compute_subtree_min_values( ... )

  // The following function computes the two smallest costs of paths
  // that connect node with index `current_index' to a leaf outside its subtree.
  // This function produces a meaningful result only if it is ran after a call
  // of the function _compute_subtree_min_values().
  template< class KernelType >
  void PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  _compute_rest_tree_min_values( Tree_type &tree, int current_index )
  {
    Node_type current_node = tree.node(current_index);

    int first_min_index = -1, second_min_index = -1;
    Number_type current_first_min(-1.0), current_second_min(-1.0);

    for( int i=0; i<current_node.marked_children.size(); i++ )
    {
      Node_type child = tree.node(current_node.marked_children[i]);
      int child_index = current_node.marked_children[i];      

      if(first_min_index == -1 || child.first_min + Number_type(child.distance) < current_first_min)
      {
        second_min_index = first_min_index;
        current_second_min = current_first_min; 
        first_min_index = child_index;
        current_first_min = child.first_min + Number_type(child.distance);
      }
      else
      if(second_min_index == -1 || child.first_min + Number_type(child.distance) < current_second_min)
      {
        second_min_index = child_index;
        current_second_min = child.first_min + Number_type(child.distance);
      }

    } // for( int i=0; i<current_node.marked_children.size(); i++ )   

    for( int i=0; i<current_node.marked_children.size(); i++ )
    {
      int child_index = current_node.marked_children[i];
      Node_type child = tree.node(current_node.marked_children[i]);

      Number_type cmin = std::min(current_first_min, current_node.rest_tree_min);

      if(current_node.rest_tree_min == Number_type(-1.0))
        cmin = current_first_min;

      if( child_index == first_min_index)
      { 
        if(current_second_min == Number_type(-1.0) || 
           (current_node.rest_tree_min < current_second_min && current_node.rest_tree_min != Number_type(-1.0)) )
          tree.node(current_node.marked_children[i]).rest_tree_min = current_node.rest_tree_min 
                                                                     + Number_type(child.distance);   
        else
          tree.node(current_node.marked_children[i]).rest_tree_min = 
                                                     current_second_min + Number_type(child.distance);   
      }
      else
        tree.node(current_node.marked_children[i]).rest_tree_min = cmin + Number_type(child.distance);

      _compute_rest_tree_min_values(tree, current_node.marked_children[i]);

    } // for( int i=0; i<current_node.marked_children.size(); i++ )

  } // _compute_rest_tree_min_values( ... )


  template< class KernelType >
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  operator()( RangeIterator rbegin, RangeIterator rend,
              int min_index, int max_index )
  {
    if(p_tree->number_of_nodes()<2)
      return Number_type(0.0);

    if(int(rend-rbegin) < 2)
      return Number_type(0.0);

    int intersection_index = p_tree->compute_intersection_node_index(min_index, max_index);

    // If the two paths coincide (intersection_index designates a leaf node)
    // then return zero distance.
    if( p_tree->node(intersection_index).children.size() == 0 )
      return Number_type(0.0);

    p_tree->node(intersection_index).mark = true;

    // Mark all the nodes which fall on a path that connects
    // a leaf node in the sample and the node indicated by intersection_index.

    p_tree->mark_Steiner_tree_of_sample(rbegin, rend);

    _compute_subtree_min_values( *p_tree, intersection_index );
    _compute_rest_tree_min_values( *p_tree, intersection_index );

    Number_type total_dist(0.0);

    // Unmark all marked nodes.

    for( RangeIterator rit = rbegin; rit != rend; rit++ )
    {
      total_dist = total_dist + p_tree->node((*rit)).rest_tree_min;

      p_tree->node(*rit).mark = false;
      p_tree->node(*rit).marked_children.clear();
      p_tree->node(*rit).marked_subtree_leaves=0;
      p_tree->node(*rit).first_min = Number_type(-1.0);
      p_tree->node(*rit).second_min = Number_type(-1.0);
      p_tree->node(*rit).rest_tree_min = Number_type(-1.0);
      Node_type v = p_tree->node(*rit);

      while( (!p_tree->is_root(v)) && p_tree->node(v.parent).mark == true )
      {
        p_tree->node(v.parent).mark = false;
        p_tree->node(v.parent).marked_children.clear();
        p_tree->node(v.parent).marked_subtree_leaves=0;
        p_tree->node(v.parent).first_min = Number_type(-1.0);
        p_tree->node(v.parent).second_min = Number_type(-1.0);
        p_tree->node(v.parent).rest_tree_min = Number_type(-1.0);
        v = p_tree->node(v.parent);

      } //while( (!p_tree->is_root(v)) && p_tree->node(v.parent).mark == true )

    } // for( RangeIterator rit = rbegin; rit != rend; rit++ )

    p_tree->clear_marked_nodes();

    return total_dist/Number_type(rend - rbegin);

  } // Mean_nearest_taxon_distance<KernelType>::operator()(...)

  template< class KernelType >
  template < class OutputIterator >
  void PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  incremental_operator_non_ultrametric( std::vector<int> &sample,
                                        std::vector<int> &sample_sizes, OutputIterator ot )
  {
    int n_l=p_tree->number_of_leaves(),
        n_n=p_tree->number_of_nodes();  

    for(int i=0; i<sample_sizes.size(); i++)
      if(sample_sizes[i] > n_l || sample_sizes[i] < 0 || sample_sizes[i] > sample.size())
      {
        std::string exception_msg;
        exception_msg += " One of the sample sizes given to the incremental operator is out of range.\n";     
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

    // If the two paths coincide (intersection_index designates a leaf node)
    // then return zero distance.
    if( p_tree->node(intersection_index).children.size() == 0 )
      *ot++ = Number_type(0.0);
    
    p_tree->node(intersection_index).mark = true;

    // Mark all the nodes which fall on a path that connects
    // a leaf node in the sample and the node indicated by intersection_index.

    typename std::vector<int>::iterator rit=sample.begin();

    for(int i=0; i<sample_sizes[ss_index]; i++)
      rit++;

    p_tree->mark_Steiner_tree_of_sample(sample.begin(), rit);
    p_tree->assign_marked_subtree_leaves(intersection_index);

    _compute_subtree_min_values( *p_tree, intersection_index );
    _compute_rest_tree_min_values( *p_tree, intersection_index );

    this->initialize_max_subtree_path_costs(intersection_index);

    Number_type total_dist(0.0);

    for( int i = 0;  i< sample_sizes[ss_index]; i++ )
      total_dist = total_dist + p_tree->node(sample[i]).rest_tree_min;

    if( p_tree->node(intersection_index).children.size() > 0 )
      *ot++ = total_dist/Number_type(sample_sizes[ss_index]);    

    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    // Continue with the rest sample sizes

    prev_sample_size = sample_sizes[ss_index];

    for(int i=ss_index+1; i<sample_sizes.size(); i++)
    {
      curr_sample_size = sample_sizes[i];

      for(int j=prev_sample_size; j<curr_sample_size; j++)
        this->update_shortest_path_costs(intersection_index,sample[j],total_dist);       
 
      *ot++ = total_dist/Number_type(curr_sample_size);
  
      prev_sample_size = curr_sample_size;

      int count = 0;

    } // for(int i=ss_index+1; i<sample_sizes.size(); i++)

    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    // Unmark all marked nodes.

    for( typename std::vector<int>::iterator rit = sample.begin(); rit != sample.end(); rit++ )
    {
      p_tree->node(*rit).mark = false;
      p_tree->node(*rit).marked_children.clear();
      p_tree->node(*rit).marked_subtree_leaves=0;
      p_tree->node(*rit).first_min = Number_type(-1.0);
      p_tree->node(*rit).second_min = Number_type(-1.0);
      p_tree->node(*rit).rest_tree_min = Number_type(-1.0);
      _max_subtree_path_costs[*rit] = Number_type(0.0);
      Node_type v = p_tree->node(*rit);

      while( (!p_tree->is_root(v)) && p_tree->node(v.parent).mark == true )
      {
        p_tree->node(v.parent).mark = false;
        p_tree->node(v.parent).marked_children.clear();
        p_tree->node(v.parent).marked_subtree_leaves=0;
        p_tree->node(v.parent).first_min = Number_type(-1.0);
        p_tree->node(v.parent).second_min = Number_type(-1.0);
        p_tree->node(v.parent).rest_tree_min = Number_type(-1.0);
        _max_subtree_path_costs[v.parent] = Number_type(0.0);
        v = p_tree->node(v.parent);

      } //while( (!p_tree->is_root(v)) && p_tree->node(v.parent).mark == true )

    } // for( std::vector<int>::iterator rit = sample.begin(); rit != sample.end(); rit++ )

    p_tree->clear_marked_nodes();

    return;

  } // Mean_nearest_taxon_distance<KernelType>::incremental_operator_non_ultrametric(...)


  template< class KernelType >
  template < class OutputIterator >
  void PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  incremental_operator_ultrametric( std::vector<int> &sample,
                                    std::vector<int> &sample_sizes, OutputIterator ot )
  {
    int n_l=p_tree->number_of_leaves(),
        n_n=p_tree->number_of_nodes();  

    for(int i=0; i<sample_sizes.size(); i++)
      if(sample_sizes[i] > n_l || sample_sizes[i] < 0 || sample_sizes[i] > sample.size())
      {
        std::string exception_msg;
        exception_msg += " One of the sample sizes given to the incremental operator is out of range.\n";     
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

    // If the two paths coincide (intersection_index designates a leaf node)
    // then return zero distance.
    if( p_tree->node(intersection_index).children.size() == 0 )
      *ot++ = Number_type(0.0);
    
    p_tree->node(intersection_index).mark = true;

    // Remember that for ultrametric trees, the value of the MNTD is
    // equal to twice the sum of weights of those edges that have exactly
    // one marked leaf (element in the sample) in their subtree.
    // Next, mark all the nodes which fall on a path that connects
    // a leaf node in the sample and the node indicated by intersection_index.
    // Then, traverse each marked edge and check if it has exactly one marked leaf 
    // in its subtree, and if so then add twice its cost in the result.

    std::vector<int>::iterator rit=sample.begin();

    for(int i=0; i<sample_sizes[ss_index]; i++)
      rit++;
    
    Number_type total_dist(0.0);
  
    p_tree->mark_Steiner_tree_of_sample(sample.begin(), rit);
    p_tree->assign_marked_subtree_leaves(intersection_index);

    if( p_tree->node(intersection_index).children.size() != 0 )
      for(int i=0; i<p_tree->number_of_marked_nodes(); i++)
        if( p_tree->node(p_tree->marked_node(i)).marked_subtree_leaves == 1 )    
          total_dist += Number_type(2.0)*p_tree->node(p_tree->marked_node(i)).distance;

    if( p_tree->node(intersection_index).children.size() > 0 )
      *ot++ = total_dist/Number_type(sample_sizes[ss_index]);

    // Continue with the rest sample sizes

    prev_sample_size = sample_sizes[ss_index];

    for(int i=ss_index+1; i<sample_sizes.size(); i++)
    {
      curr_sample_size = sample_sizes[i];

      for(int j=prev_sample_size; j<curr_sample_size; j++)
        total_dist += this->update_total_cost_ultrametric(intersection_index,sample[j]); // intersection_index
                                                                                          // might get 
                                                                                          // updated here       

      *ot++ = total_dist/Number_type(curr_sample_size);
  
      prev_sample_size = curr_sample_size;

    } // for(int i=ss_index+1; i<sample_sizes.size(); i++)


    // Unmark all marked nodes.

    p_tree->unmark_Steiner_tree_of_sample(sample.begin(), sample.end());

    return;

  } // Mean_nearest_taxon_distance<KernelType>::incremental_operator_ultrametric(...)


  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>:: 
  update_total_cost_ultrametric(int &intersection_index, int new_node_index)
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

        previous_index = current_index;
        current_index = p_tree->node(current_index).parent;
        previous_node = p_tree->node(previous_index);

      }while(previous_index != intersection_index);

    } // if(old_intersection_index != intersection_index)

    Number_type total_added_dist(0.0);

    current_index = p_tree->node(new_node_index).parent;
    previous_index = new_node_index; 
    previous_node = p_tree->node(new_node_index);

    p_tree->node(new_node_index).marked_subtree_leaves=1;

    total_added_dist += Number_type(2.0)*previous_node.distance;

    p_tree->insert_marked_node(new_node_index);
    p_tree->node(new_node_index).mark = true;

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
  
      p_tree->node(current_index).marked_subtree_leaves++;

      if(p_tree->node(current_index).marked_subtree_leaves==2)
        total_added_dist -= Number_type(2.0)*p_tree->node(current_index).distance;
      else if(p_tree->node(current_index).marked_subtree_leaves==1)
        total_added_dist += Number_type(2.0)*p_tree->node(current_index).distance;
      else
        break;

      previous_index = current_index;
      current_index = p_tree->node(current_index).parent;
      previous_node = p_tree->node(previous_index);

    }while(previous_index != intersection_index);

    return total_added_dist;

  } // update_total_cost_ultrametric(...)

  template< class KernelType >
  void PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::  
  initialize_max_subtree_path_costs(int index)
  {
    if(p_tree->node(index).number_of_children() == 0)
    { 
       this->_max_subtree_path_costs[index] = p_tree->node(index).rest_tree_min;
       return; 
    }

    Number_type max(-1.0);

    for(int i=0; i<p_tree->node(index).number_of_marked_children(); i++)
    {
      int child = p_tree->node(index).marked_children[i];

      initialize_max_subtree_path_costs(child); 

      if( max < Number_type(0.0) || this->_max_subtree_path_costs[child] > max)
        max = this->_max_subtree_path_costs[child];

    } // for(int i=0; i<p_tree->node(index).number_of_marked_children(); i++)

    this->_max_subtree_path_costs[index] = max;

  } // initialize_max_subtree_path_costs(int index)

  template< class KernelType >
  void PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  update_shortest_path_costs(int &intersection_index,int new_node_index, Number_type &total_dist)
  {
     if(intersection_index == new_node_index)
      return;

    int old_intersection_index = intersection_index;
    int current_index, previous_index;
    Node_type previous_node;

    p_tree->update_marked_Steiner_tree(intersection_index, new_node_index);

    if(old_intersection_index != intersection_index)
    {
      int current_index = p_tree->node(old_intersection_index).parent,
          previous_index = old_intersection_index;

      previous_node = p_tree->node(previous_index);

      do
      {
        this->_max_subtree_path_costs[current_index]= this->_max_subtree_path_costs[previous_index] + 
                                                       previous_node.distance;

        p_tree->node(current_index).marked_subtree_leaves = previous_node.marked_subtree_leaves;
        p_tree->node(current_index).first_min = previous_node.first_min + previous_node.distance;
        p_tree->insert_marked_node(current_index);

        previous_index = current_index;
        current_index = p_tree->node(current_index).parent;
        previous_node = p_tree->node(previous_index);

      }while(previous_index != intersection_index);

    } // if(old_intersection_index != intersection_index)


    Number_type path_distance=0.0;

    current_index = p_tree->node(new_node_index).parent;
    previous_index = new_node_index; 
    previous_node = p_tree->node(new_node_index);

    p_tree->node(new_node_index).marked_subtree_leaves=1;
    p_tree->node(new_node_index).first_min = Number_type(0.0);
    p_tree->node(new_node_index).mark = true;
    p_tree->insert_marked_node(new_node_index);

    int min_index=-1; 
    Number_type min_distance(-1.0);

    std::vector<int> NN_leaves;

    // We have to insert to the marked_nodes vector the new nodes that appear 
    // in the tree and which are exclusive ancestors of the newly added leaf.
    // These fall inside the maximal chain of nodes that have degree exactly 
    // two and which connects the new leaf with rest of the tree.
    // Thus, the last of the nodes in this chain is connected with an edge
    // to a node with degree > 2. As soon as we find this condition, we stop
    // inserting nodes to the marked nodes vector.

    bool found_edge_to_existing_tree=false;

    Number_type distance_difference(0.0); 

    do
    {
      // The next condition checks if the current index points to a 
      // node that has been freshly added to the marked subtree,
      // and it is different from the root of this subtree.

      if(found_edge_to_existing_tree==false && p_tree->node(current_index).number_of_marked_children()==1)
        p_tree->insert_marked_node(current_index);
      else if(found_edge_to_existing_tree==false)
        found_edge_to_existing_tree=true;

      p_tree->node(current_index).marked_subtree_leaves++;
      path_distance += previous_node.distance;

      if( p_tree->node(current_index).number_of_marked_children() > 1 &&
          ( p_tree->node(current_index).first_min+path_distance < min_distance || 
            min_distance < Number_type(0.0) ) )
      {
        min_index = current_index;
        min_distance = p_tree->node(current_index).first_min+path_distance;
      }

      if(p_tree->node(current_index).first_min > path_distance || p_tree->node(current_index).first_min  < Number_type(0.0))
        p_tree->node(current_index).first_min = path_distance;     
      
      // The newly added leaf may be the new nearest neighbour 
      // for marked leaf nodes that lie in this subtree.
      // Next we check this out.

      for(int i=0; i<p_tree->node(current_index).number_of_marked_children(); i++)
      {
        int curr_child = p_tree->node(current_index).marked_children[i];

        if( curr_child != previous_index &&
            path_distance + p_tree->node(curr_child).distance < this->_max_subtree_path_costs[curr_child])
          this->find_new_nearest_neighbours(path_distance + p_tree->node(curr_child).distance,
                                            curr_child, std::back_inserter(NN_leaves), distance_difference);
      }

      // Done with this node, go one node up.

      previous_index = current_index;
      current_index = p_tree->node(current_index).parent;
      previous_node = p_tree->node(previous_index);

    }while(previous_index != intersection_index);

    p_tree->node(new_node_index).rest_tree_min = min_distance;
    this->_max_subtree_path_costs[new_node_index] = min_distance;
    total_dist += p_tree->node(new_node_index).rest_tree_min;
    total_dist += distance_difference;

    NN_leaves.push_back(new_node_index);
 
    for(int i=0; i<NN_leaves.size(); i++)
      this->update_max_subtree_path_costs(NN_leaves[i]);    
  
  } // update_shortest_path_costs(...)

  template< class KernelType >
  template< class OutputIterator >  
  void PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::  
  find_new_nearest_neighbours(Number_type dist, int index, OutputIterator ot, 
                              Number_type &distance_difference)
  {
    if(p_tree->node(index).number_of_children() == 0)
    {
      distance_difference += dist-p_tree->node(index).rest_tree_min; 
      p_tree->node(index).rest_tree_min = dist;
      this->_max_subtree_path_costs[index] = dist;
      *ot++ = index;
      return;
    }

    for(int i=0; i<p_tree->node(index).number_of_marked_children(); i++)
    {
      int curr_child = p_tree->node(index).marked_children[i];

      if( dist + p_tree->node(curr_child).distance < this->_max_subtree_path_costs[curr_child])
        this->find_new_nearest_neighbours(dist + p_tree->node(curr_child).distance, 
                                          curr_child, ot, distance_difference);
    }

  } // find_new_nearest_neighbours(...)

  template< class KernelType >
  void PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::  
  update_max_subtree_path_costs(int index)
  {
     if(p_tree->node(index).mark == false)
       return;

     if(p_tree->node(index).number_of_children() == 0 && index != p_tree->root_index())
       update_max_subtree_path_costs(p_tree->node(index).parent);
     else
     {
       Number_type max(-1.0);

       for(int i=0; i<p_tree->node(index).number_of_marked_children(); i++)
       { 
         if( max < Number_type(0.0) || 
             this->_max_subtree_path_costs[p_tree->node(index).marked_children[i]] > max)
           max = this->_max_subtree_path_costs[p_tree->node(index).marked_children[i]];
       }
 
       this->_max_subtree_path_costs[index] = max;

       if(index != p_tree->root_index())
         update_max_subtree_path_costs(p_tree->node(index).parent); 

     } // else of if(p_tree->node(index).number_of_children() == 0 && ... )

  } // update_max_subtree_path_costs(int index)

  template< class KernelType >
  template< class OutputIterator >
  void PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  _compute_subtree_sums( int index, Number_type& sum_of_products, OutputIterator ot,
                         Number_type &sum_subtree, Number_type &sum_subtract )
  {
    Node_type v = p_tree->node(index);

    Number_type mhyperg = this->hypergeom_minus_one(_number_of_leaves-v.all_subtree_leaves);

    std::vector<int> subtree_leaves;

    for( int i=0; i<v.children.size(); i++ )
    {
      Number_type sum_of_pr(0.0);
      std::vector< std::pair<Number_type,int> > c_subtree_leaves;

      _compute_subtree_sums( v.children[i], sum_of_pr,
                             std::back_inserter(c_subtree_leaves),
                             sum_subtree, sum_subtract );

      sum_subtree  += Number_type(v.distance)*sum_of_pr*mhyperg;

      for( int j=0; j< c_subtree_leaves.size(); j++ )
      {
        sum_subtract += Number_type(v.distance)*c_subtree_leaves[j].first*
                        Number_type(v.all_subtree_leaves)*c_subtree_leaves[j].second*
                        hypergeom_minus_two(_number_of_leaves-v.all_subtree_leaves-c_subtree_leaves[j].second);

        *ot++ = c_subtree_leaves[j];
      }

      sum_of_products += sum_of_pr;
    }

    sum_subtree  += Number_type(v.distance)*Number_type(v.distance)*Number_type(v.all_subtree_leaves)*mhyperg;

    sum_subtract += Number_type(v.distance)*Number_type(v.distance)*
                    Number_type(v.all_subtree_leaves)*Number_type(v.all_subtree_leaves)*
                    two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT);

    sum_of_products += Number_type(v.distance)*Number_type(v.all_subtree_leaves);
    *ot++ = std::make_pair(Number_type(v.distance), v.all_subtree_leaves);

  } // _compute_subtree_sums( int index, ... )


  // Computes all together the probability values f(x) = \binom{x}{r}/\binom{s}{r}

  template< class KernelType >
  void PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  compute_all_hypergeometric_probabilities( int sample_size, int number_of_leaves)
  {
    _sample_size = sample_size;
    _number_of_leaves    = number_of_leaves;

    if( !_hypergeom.empty() )
      _hypergeom.clear();

    std::vector<Number_type> tempgeom;
    tempgeom.push_back(Number_type(1.0));

    for( int i= _number_of_leaves-1; i>=_sample_size; i-- )
    {
      Number_type x(i+1);
      tempgeom.push_back( tempgeom.back()/ Number_type(x/(x-Number_type(_sample_size)))  );
    }

    for( int i= tempgeom.size()-1; i>=0; i-- )
      _hypergeom.push_back( tempgeom[i] );
	  
  } // compute_all_hypergeometric_probabilities(...)


  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  compute_expectation_uniform_distribution( int sample_size )
  {
    p_tree->assign_all_subtree_leaves(p_tree->root_index());
    compute_all_hypergeometric_probabilities( sample_size, p_tree->number_of_leaves());

    Number_type expectation(0.0);

    for( int i=0; i<p_tree->number_of_nodes()-1; i++ ) // We exclude the edge of the root node.
        expectation += Number_type(p_tree->node(i).distance)*Number_type(p_tree->node(i).all_subtree_leaves)*
                       hypergeom_minus_one(p_tree->number_of_leaves()-p_tree->node(i).all_subtree_leaves);

    return expectation*Number_type(2.0)/Number_type(sample_size);
	
  } // compute_expectation_uniform_distribution( int sample_size )

  template< class KernelType >
  typename KernelType::Number_type 
  PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::compute_expectation( int sample_size )
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

    if(!p_tree->is_ultrametric())
    {
      std::string exception_msg;
      exception_msg += " Request to compute MNTD expectation on a non-ultrametric tree.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if( sample_size <= 1)
      return Number_type(0.0);

    if(this->probability_distribution() == Kernel::UNIFORM_FIXED_SIZE)
      return compute_expectation_uniform_distribution(sample_size);
    else if(this->probability_distribution() == Kernel::POISSON_BINOMIAL_FIXED_SIZE)
    {
      if(sample_size > _CPoisson_exps.size()-1 || _CPoisson_exps.size() == 0)
      {
        _CPoisson_exps.clear();
        _CPoisson_vars.clear();

        Poisson_moments_functor().compute_expectations_and_variances( *p_tree, sample_size, 
                                                                      std::back_inserter(_CPoisson_exps), 
                                                                      std::back_inserter(_CPoisson_vars) );
      }

      return _CPoisson_exps[sample_size];

    } // else if(this->probability_distribution() == Kernel::POISSON_BINOMIAL_FIXED_SIZE)
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
    {
      // (this->probability_distribution() == Kernel::POISSON_BINOMIAL);

      return Number_type(-1.0);
    }

  } // compute_expectation( int sample_size )
  
  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  compute_variance_uniform_distribution( int sample_size, Number_type expect )
  {

    if( sample_size <= 1 || sample_size == p_tree->number_of_leaves() )
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(p_tree->root_index());

    Number_type exp;

    if( expect != Number_type(-1.0) )
      exp = expect;
    else
      exp = compute_expectation( sample_size );

    Number_type sum_subtree(0.0), sum_subtract(0.0), sum_third_case(0.0),
                sum_self(0.0), sum_self_third_case(0.0), sum_same_class_third_case(0.0), total_sum(0.0);

    _compute_subtree_sums(sum_subtree, sum_subtract);

    std::map< int, Number_type > edge_classes;

    for( int i=0; i<p_tree->number_of_nodes()-1; i++ )
    {
      Node_type v = p_tree->node(i);

      if( edge_classes.find(v.all_subtree_leaves) == edge_classes.end() )
        edge_classes[v.all_subtree_leaves] = Number_type(v.distance)*Number_type(v.all_subtree_leaves);
      else
        edge_classes[v.all_subtree_leaves] += Number_type(v.distance)*Number_type(v.all_subtree_leaves);

      sum_self += Number_type(v.distance)*Number_type(v.distance)*Number_type(v.all_subtree_leaves)*
                 two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR);

      sum_self_third_case += Number_type(v.distance)*Number_type(v.distance)*
                             Number_type(v.all_subtree_leaves)*Number_type(v.all_subtree_leaves)*
                             two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT);
    }

    typename std::map< int, Number_type >::iterator cit1, cit2;

    for(cit1 = edge_classes.begin(); cit1 != edge_classes.end(); cit1++ )
    {
      sum_same_class_third_case += cit1->second*cit1->second*
                                   two_edge_pr(cit1->first, cit1->first, Kernel::INDEPENDENT);



      for(cit2 = edge_classes.begin(); cit2 != cit1; cit2++ )
        sum_third_case += cit1->second*cit2->second*
                          two_edge_pr(cit1->first, cit2->first, Kernel::INDEPENDENT);
    }

    sum_same_class_third_case = ((sum_same_class_third_case - sum_self_third_case)/Number_type(2.0))
                                  + sum_self_third_case;

    sum_third_case += sum_same_class_third_case;

    total_sum = (Number_type(2.0)*(sum_third_case - sum_subtract + sum_subtree)) - sum_self;
    total_sum = Number_type(4.0)*total_sum/(Number_type(_sample_size)*Number_type(_sample_size));

    return total_sum-(exp*exp);

  } //compute_variance_uniform_distribution( ... )

  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  compute_variance_uniform_distribution_slow( int sample_size, Number_type expect )
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

    if(!p_tree->is_ultrametric())
    {
      std::string exception_msg;
      exception_msg += " Request to compute MNTD variance on a non-ultrametric tree.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if( sample_size <= 1 )
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(*p_tree, p_tree->root_index());

    Number_type exp;

    if( expect != Number_type(-1.0) )
      exp = expect;
    else
      exp = compute_expectation( sample_size );

    if(p_tree->subtree_edges_size() == 0)
      p_tree->compute_subtree_edges(p_tree->root_index());

    Number_type sum(0.0);

    for( int i=0; i<p_tree->number_of_nodes()-1; i++ )
    {
      Node_type u = p_tree->node(i);

      for( int j=0; j<p_tree->number_of_nodes()-1; j++ )
      {
        Node_type v = p_tree->node(j);

        if( j < i && j >= i - p_tree->subtree_edges(i) )
          sum += Number_type(u.distance)*Number_type(v.distance)*
                 Number_type(v.all_subtree_leaves)*two_edge_pr(u.all_subtree_leaves, 
                                                               v.all_subtree_leaves, Kernel::OFFSPRING);
        else if( i <= j && i >= j - p_tree->subtree_edges(j) )
          sum += Number_type(u.distance)*Number_type(v.distance)*
                 Number_type(u.all_subtree_leaves)*two_edge_pr(u.all_subtree_leaves, 
                                                               v.all_subtree_leaves, Kernel::ANCESTOR);
        else
          sum += Number_type(u.distance)*Number_type(v.distance)*Number_type(u.all_subtree_leaves)*
                 Number_type(v.all_subtree_leaves)*two_edge_pr(u.all_subtree_leaves, 
                                                               v.all_subtree_leaves, Kernel::INDEPENDENT);

      } // for( int j=0; j<p_tree->number_of_nodes()-1; j++ )

    } // for( int i=0; i<p_tree->number_of_nodes()-1; i++ )


    return (Number_type(4.0)*sum/((sample_size)*(sample_size)))-(exp*exp);

  } // compute_variance_slow( ... )

  template< class KernelType >
  typename KernelType::Number_type 
  PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  compute_variance( int sample_size, Number_type expect)
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

    if(!p_tree->is_ultrametric())
    {
      std::string exception_msg;
      exception_msg += " Request to compute MNTD variance on a non-ultrametric tree.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(sample_size <= 1)
      return Number_type(0.0);

    Number_type variance;   

    if(this->probability_distribution() == Kernel::UNIFORM_FIXED_SIZE)
      variance = compute_variance_uniform_distribution(sample_size);
    else if(this->probability_distribution() == Kernel::POISSON_BINOMIAL_FIXED_SIZE)
    {
      if(sample_size > _CPoisson_vars.size()-1 || _CPoisson_vars.size() == 0)
      {
        _CPoisson_exps.clear();
        _CPoisson_vars.clear();

        Poisson_moments_functor().compute_expectations_and_variances( *p_tree, sample_size, 
                                                                      std::back_inserter(_CPoisson_exps), 
                                                                      std::back_inserter(_CPoisson_vars) );
      }

      variance = _CPoisson_vars[sample_size];
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
    {
      // (this->probability_distribution() == Kernel::POISSON_BINOMIAL);

      return Number_type(-1.0);
    }

    return variance;

  } // compute_variance(int sample_size)

  template< class KernelType >
  typename KernelType::Number_type 
  PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  compute_deviation( int sample_size, Number_type expect)
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

    if(!p_tree->is_ultrametric())
      return Number_type(-1.0);

    Number_type variance;   

    if(this->probability_distribution() == Kernel::UNIFORM_FIXED_SIZE)
      variance = compute_variance_uniform_distribution(sample_size, expect);
    else if(this->probability_distribution() == Kernel::POISSON_BINOMIAL_FIXED_SIZE)
    {
      if(sample_size > _CPoisson_vars.size()-1 || _CPoisson_vars.size() == 0)
      {
        _CPoisson_exps.clear();
        _CPoisson_vars.clear();

        Poisson_moments_functor().compute_expectations_and_variances( *p_tree, sample_size, 
                                                                      std::back_inserter(_CPoisson_exps), 
                                                                      std::back_inserter(_CPoisson_vars) );
      }

      variance = _CPoisson_vars[sample_size];
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
    {
      // (this->probability_distribution() == Kernel::POISSON_BINOMIAL);

      return Number_type(-1.0);
    }
 
    if( variance < Number_type(0.0) ) 
      return Number_type(0.0);

    return Square_root()(variance); 

  } // compute_deviation(...)

#endif //MEAN_NEAREST_TAXON_DISTANCE_IMPL_H
