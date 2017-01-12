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

#ifndef CORE_ANCESTOR_COST_IMPL_H
#define CORE_ANCESTOR_COST_IMPL_H

  template< class KernelType >
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  operator()( RangeIterator rbegin, RangeIterator rend,
              int min_index, int max_index)
  {
    if(p_tree->number_of_nodes()<2)
      return Number_type(0.0);

    int number_of_input_leaves = rend-rbegin;
    Number_type temp_rchi = Ceiling()(this->chi()*Number_type(number_of_input_leaves)); 
    int rchi = int(To_double()(temp_rchi));

    if(rchi == 0)
      return Number_type(0.0);

    // Mark all the nodes which fall on a path that connects
    // a leaf node in the input sample and the root node.

    p_tree->mark_Steiner_tree_of_sample(rbegin, rend);
    p_tree->assign_marked_subtree_leaves(p_tree->root_index());

    // Starting from the root, traverse the tree to find
    // the core ancestor of the input leaf sample.

    Number_type total_dist = Number_type(0.0);
    Node_type v = p_tree->root();
    bool found = false;

    int last_heir = p_tree->root_index();

    do
    {
      if(v.number_of_marked_children() > 0)
      {
        int heir=-1;

        for(int i=0; i<v.number_of_marked_children(); i++)
          if( p_tree->node(v.marked_children[i]).marked_subtree_leaves >= rchi )
          {
            heir = v.marked_children[i];
	    break;
          }

	if(heir == -1)
	  found = true;
	else
	{
	  total_dist += Number_type(p_tree->node(heir).distance);
	  v = p_tree->node(heir);
	  last_heir = heir;
	}

      } // if(v.number_of_children() > 0)

    }while(found == false && v.number_of_children() > 0);

    p_tree->unmark_Steiner_tree_of_sample(rbegin, rend);

    return total_dist;

  } // Core_ancestor_cost<KernelType>::operator()(...)


  template< class KernelType >
  template < class OutputIterator >
  void PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
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

    while(ss_index < sample_sizes.size() && sample_sizes[ss_index] <1)
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

    Number_type temp_rchi = Ceiling()(this->chi()*Number_type(sample_sizes[ss_index])); 
    int rchi = int(To_double()(temp_rchi));
    
    // Mark all the nodes which fall on a path that connects
    // a leaf node in the sample and the root.
    // For each edge that we traverse, add the corresponding weight to the solution.

    std::vector<int>::iterator rit=sample.begin();

    for(int i=0; i<sample_sizes[ss_index]; i++)
      rit++;
      
    p_tree->mark_Steiner_tree_of_sample(sample.begin(), rit);
    p_tree->assign_marked_subtree_leaves(p_tree->root_index());

    // Starting from the root, traverse the tree to find
    // the core ancestor of the input leaf sample.

    Number_type total_dist = Number_type(0.0);
    Node_type v = p_tree->root();
    bool found = false;
    int last_heir = p_tree->root_index();

    do
    {
      if(v.number_of_marked_children() > 0)
      {
        int heir=-1;

        for(int i=0; i<v.number_of_marked_children(); i++)
          if( p_tree->node(v.marked_children[i]).marked_subtree_leaves >= rchi )
          {
            heir = v.marked_children[i];
	    break;
          }

	if(heir == -1)
	  found = true;
	else
	{
	  total_dist += Number_type(p_tree->node(heir).distance);
	  v = p_tree->node(heir);
	  last_heir = heir;
	}

      } // if(v.number_of_children() > 0)

    }while(found == false && v.number_of_children() > 0);

    *ot++ = total_dist;

    // Continue with the rest sample sizes

    prev_sample_size = sample_sizes[ss_index];

    for(int i=ss_index+1; i<sample_sizes.size(); i++)
    {
      curr_sample_size = sample_sizes[i];

      temp_rchi = Ceiling()(this->chi()*Number_type(curr_sample_size)); 
      rchi = int(To_double()(temp_rchi));

      for(int j=prev_sample_size; j<curr_sample_size; j++)
         total_dist = this->update_marked_subtree(sample[j],rchi);

      if(rchi == 0)
	*ot++ = Number_type(0.0);
      else
        *ot++ = total_dist;
  
      prev_sample_size = curr_sample_size;

    } // for(int i=ss_index+1; i<sample_sizes.size(); i++)

    // Unmark all marked nodes.

    p_tree->unmark_Steiner_tree_of_sample(sample.begin(), sample.end());

    return;

  } // incremental_operator(...)

  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  update_marked_subtree(int new_node_index, int rchi)
  {
    int current_index, previous_index;
    Node_type previous_node;

    Number_type path_dist=0.0;

    current_index = p_tree->node(new_node_index).parent;
    previous_index = new_node_index; 
    previous_node = p_tree->node(new_node_index);

    p_tree->node(new_node_index).marked_subtree_leaves=1;
    p_tree->node(new_node_index).mark = true;

    // We have to update the marked_children member vector for all nodes
    // in the tree which are exclusive ancestors of the newly added leaf.
    // These fall inside the maximal chain of nodes that have degree exactly 
    // two and which connects the new leaf with rest of the tree.
    // Thus, the last of the nodes in this chain is connected with an edge
    // to a node with degree > 2. 

    bool found_edge_to_existing_tree=false;

    if(current_index != -1)
      do
      { 
        p_tree->node(current_index).marked_subtree_leaves++;
        p_tree->node(current_index).mark = true;


        if(found_edge_to_existing_tree==false && p_tree->node(current_index).number_of_marked_children()==0)
          p_tree->node(current_index).marked_children.push_back(previous_index);
        else if(found_edge_to_existing_tree==false)
        {
          p_tree->node(current_index).marked_children.push_back(previous_index);
          found_edge_to_existing_tree=true;
        }

        previous_index = current_index;
        current_index = p_tree->node(current_index).parent;
        previous_node = p_tree->node(previous_index);

      }while(previous_index != p_tree->root_index());

    Number_type total_dist = Number_type(0.0);
    int last_heir = p_tree->root_index();
    Node_type v = p_tree->root();
    bool found = false;

    do
    {
      if(v.number_of_marked_children() > 0)
      {
        int heir=-1;

        for(int i=0; i<v.number_of_marked_children(); i++)
          if( p_tree->node(v.marked_children[i]).marked_subtree_leaves >= rchi )
          {
            heir = v.marked_children[i];
	    break;
          }

	if(heir == -1)
	  found = true;
	else
	{
	  total_dist += Number_type(p_tree->node(heir).distance);
	  v = p_tree->node(heir);
	  last_heir = heir;
	}

      } // if(v.number_of_children() > 0)

    }while(found == false && v.number_of_children() > 0);

    return total_dist;

  } // update_marked_subtree(...)

  template< class KernelType >
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  slow_operator( RangeIterator rbegin, RangeIterator rend,
                 int min_index, int max_index)
  {
    int number_of_input_leaves = rend-rbegin;
    Number_type temp_rchi = Ceiling()(this->chi()*Number_type(number_of_input_leaves));
	int rchi = int(To_double(temp_rchi));

    if(rchi == 0)
	  return Number_type(0.0);

    for(int i=0; i<p_tree->number_of_nodes(); i++)
    {
      p_tree->node(i).marked_subtree_leaves = 0;
      p_tree->node(i).mark= false;
    }

    for(RangeIterator rit=rbegin; rit != rend; rit++)
    {
      p_tree->node(*rit).mark=true;

      Node_type u = p_tree->node(*rit);

      while(u.parent >=0)
      {
        p_tree->node(u.parent).mark=true;
        u = p_tree->node(u.parent);
      }
    }

    p_tree->assign_marked_subtree_leaves(p_tree->root_index());

    int largest_index=-1;

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      if(p_tree->node(i).marked_subtree_leaves >= rchi && i>largest_index  )
      {
        int new_index= -1;

        for(int j=0; j<p_tree->node(i).number_of_children(); j++)
        {
          int child = p_tree->node(i).children[j];

          if(p_tree->node(child).marked_subtree_leaves >= rchi)
            new_index = child;
        }

        if(new_index==-1)
          largest_index=i;

      } // if(p_tree->node(i).marked_subtree_leaves >= ... )

    Node_type v= p_tree->node(largest_index);
    Number_type dist(0.0);

    while(v.parent >=0)
    {
      dist = dist + Number_type(v.distance);
      v = p_tree->node(v.parent);
    }

    for(int i=0; i<p_tree->number_of_nodes(); i++)
    {
      p_tree->node(i).marked_subtree_leaves = 0;
      p_tree->node(i).mark= false;
    }

    return dist;

  } //slow_operator(...)

  // Computes all together the probability values f(x) = \binom{x}{r}/\binom{s}{r}

  template< class KernelType >
  void PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_all_hypergeometric_probabilities_a( int sample_size, int number_of_leaves)
  {
    _sample_size = sample_size;
    _number_of_leaves = number_of_leaves;

    if( !_hypergeom_a.empty() )
      _hypergeom_a.clear();

    std::vector<Protected_number_type> tempgeom;
    tempgeom.push_back(Number_type(1.0));

    for( int i= _number_of_leaves-1; i>=_sample_size; i-- )
    {
      Protected_number_type x(Number_type(i+1));

      tempgeom.push_back( tempgeom.back()*Protected_number_type(Number_type(i+1-_sample_size))/x  );
    }

    for( int i= tempgeom.size()-1; i>=0; i-- )
      _hypergeom_a.push_back( tempgeom[i] );

  } //compute_all_hypergeometric_probabilities_a(...)


  // Computes all together the probability values f(x) = \binom{x}{s-r}/\binom{s}{r}

  template< class KernelType >
  void PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_all_hypergeometric_probabilities_b( int sample_size, int number_of_leaves)
  {
    _sample_size = sample_size;
    _number_of_leaves = number_of_leaves;

    if( !_hypergeom_b.empty() )
      _hypergeom_b.clear();

    std::vector<Protected_number_type> tempgeom;
    tempgeom.push_back(Protected_number_type(Number_type(1.0)));

    for( int i=_number_of_leaves-1; i>=_number_of_leaves-_sample_size; i-- )
    {
      Protected_number_type x(Number_type(i+1));

      tempgeom.push_back( tempgeom.back()*
                          Protected_number_type(Number_type(_sample_size+i+1-_number_of_leaves))/x );
    }

    for( int i= tempgeom.size()-1; i>=0; i-- )
      _hypergeom_b.push_back( tempgeom[i] );

  } //compute_all_hypergeometric_probabilities_b(...)


  template< class KernelType >
  typename KernelType::Numeric_traits::Protected_number_type
  PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_node_probability(int number_of_subtree_leaves, int sample_size, bool silent)
  {
    int number_of_all_leaves = p_tree->number_of_leaves();
    
    Number_type temp_rchi = Ceiling()(_chi*Number_type(sample_size));
    int leaf_fraction = int(To_double()(temp_rchi));

    if(leaf_fraction < 1)
      return Protected_number_type(Number_type(0.0));

    if(std::max(1,sample_size+number_of_subtree_leaves-number_of_all_leaves) > 
                  std::min(number_of_subtree_leaves,sample_size))
      return Protected_number_type(Number_type(0.0));

    if(leaf_fraction > number_of_subtree_leaves)
      return Protected_number_type(Number_type(0.0));

    if(number_of_subtree_leaves == number_of_all_leaves)
      return Protected_number_type(Number_type(1.0));

    Protected_number_type basic_hypergeom, extra_hypergeom(Number_type(0.0)),
		          probability_product(Number_type(1.0)),
		          cumulative_probability(Number_type(0.0));

    int i;

    if( number_of_all_leaves-sample_size-number_of_subtree_leaves >= 0 )
    {
      basic_hypergeom = hypergeom_a(number_of_all_leaves - number_of_subtree_leaves);
      i=1;
    }
    else
    {
      basic_hypergeom = hypergeom_b(number_of_subtree_leaves);
      i=sample_size+number_of_subtree_leaves-number_of_all_leaves+1;
    }

    if(i-1>=leaf_fraction)
      extra_hypergeom = basic_hypergeom;

    probability_product = basic_hypergeom;

    for( /*We settled i already*/ ; i<=std::min(number_of_subtree_leaves,sample_size); i++ )
    {
      Number_type f1(sample_size-i+1), f2(number_of_subtree_leaves-i+1),
	          f3(std::max(number_of_all_leaves-number_of_subtree_leaves-sample_size+i,1)), 
                  f4(i); 

      Protected_number_type numerator, denominator;

      numerator = Protected_number_type(Protected_number_type(f1)*Protected_number_type(f2));
      denominator = Protected_number_type(Protected_number_type(f3)*Protected_number_type(f4));
      probability_product = probability_product*(numerator/denominator);

      if(i>=leaf_fraction)
        cumulative_probability+=probability_product;
    }

    //return basic_hypergeom*(cumulative_probability+extra_hypergeom);

    return cumulative_probability+extra_hypergeom;

  } // compute_node_probability(...)


  template< class KernelType >
  template<class OutputIterator>
  void PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_all_root_path_costs(OutputIterator ot)
  {
    std::vector<Number_type> path_costs;
    path_costs.assign( p_tree->number_of_nodes() , Number_type(0.0) );

    std::queue< std::pair<int,Number_type> > node_cost_pairs;
	node_cost_pairs.push( std::make_pair(p_tree->root_index(), Number_type(0.0)) );

    while(!node_cost_pairs.empty())
    {
      int index = node_cost_pairs.front().first;
	  Number_type cost = node_cost_pairs.front().second;

      node_cost_pairs.pop();
      path_costs[index] = cost;

      Node_type v = p_tree->node(index);

      for(int i=0; i<v.number_of_children(); i++)
      {
	int child_index = v.children[i];
        Number_type new_cost = cost + Number_type(p_tree->node(child_index).distance);
        node_cost_pairs.push(std::make_pair(child_index,new_cost));
      }

    } // while(!node_indices.empty())

    for(int i=0; i<path_costs.size(); i++)
      *ot++ = Protected_number_type(path_costs[i]);

  } // compute_all_root_path_costs(OutputIterator ot)


  template< class KernelType >
  template<class OutputIterator>
  void PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_first_k_raw_moments_protected(int k, int sample_size, OutputIterator ot)
  {
    if(sample_size < 0 || sample_size > p_tree->number_of_leaves())
    {
      std::string exception_msg;
      exception_msg += " Request to compute moments with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(p_tree->number_of_leaves()<=1 || sample_size == 0)
    {
      for(int i=0; i<k; i++)
        *ot++ = Protected_number_type(Number_type(0.0));

      return;
    }

    std::vector<Protected_number_type> root_path_costs, subtree_size_probabilities,
	                               node_probabilities, costs_kth_exponent;

    p_tree->assign_all_subtree_leaves(p_tree->root_index());

    this->compute_all_hypergeometric_probabilities_a(sample_size, p_tree->number_of_leaves());
    this->compute_all_hypergeometric_probabilities_b(sample_size, p_tree->number_of_leaves());
    this->compute_all_root_path_costs(std::back_inserter(root_path_costs));
    subtree_size_probabilities.assign(p_tree->number_of_leaves(), Protected_number_type(Number_type(0.0)));
    node_probabilities.assign(p_tree->number_of_nodes(), Protected_number_type(Number_type(0.0)));
    costs_kth_exponent.assign(p_tree->number_of_nodes(), Protected_number_type(Number_type(1.0)));

    std::set<int> unique_subtree_sizes;

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      unique_subtree_sizes.insert(p_tree->node(i).all_subtree_leaves);

    typename std::set<int>::iterator sit;

    for( sit=unique_subtree_sizes.begin(); sit!=unique_subtree_sizes.end(); sit++)
    {
      subtree_size_probabilities[(*sit)-1] = this->compute_node_probability(*sit,sample_size);

      ////////////////////////////////////////////////////////////////////////////
      //if(this->compute_node_probability(*sit,sample_size) > Number_type(1.0)) //
      //  this->compute_node_probability(*sit,sample_size,false);               // 
      ////////////////////////////////////////////////////////////////////////////
    }

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      node_probabilities[i] = subtree_size_probabilities[p_tree->node(i).all_subtree_leaves-1];

    std::vector<Protected_number_type> final_probabilities = node_probabilities;

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      for( int j=0; j<p_tree->node(i).number_of_children(); j++ )
        final_probabilities[i] = final_probabilities[i] - node_probabilities[p_tree->node(i).children[j]];

    for(int i=0; i<k; i++)
    {
      Protected_number_type ith_raw_moment(0.0);

      for(int j=0; j < p_tree->number_of_nodes(); j++)
      {
        costs_kth_exponent[j] = costs_kth_exponent[j]*root_path_costs[j];
        ith_raw_moment += costs_kth_exponent[j]*final_probabilities[j];
      }

      *ot++ = ith_raw_moment;

    } // for(int i=0; i<k; i++)

  } //compute_first_k_raw_moments_protected(int k, OutputIterator ot)



  template< class KernelType >
  template<class OutputIterator>
  void PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_first_k_centralised_moments(int k, int sample_size, OutputIterator ot)
  {
    if(sample_size < 0 || sample_size > p_tree->number_of_leaves())
    {
      std::string exception_msg;
      exception_msg += " Request to compute moments with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(p_tree->number_of_leaves()<=1 || sample_size == 0)
    {
      for(int i=0; i<k; i++)
        *ot++ = Protected_number_type(Number_type(0.0)).to_number_type();

      return;
    }

    std::vector<Protected_number_type> root_path_costs, subtree_size_probabilities,
	                               node_probabilities, costs_kth_exponent;

    p_tree->assign_all_subtree_leaves(p_tree->root_index());

    this->compute_all_hypergeometric_probabilities_a(sample_size, p_tree->number_of_leaves());
    this->compute_all_hypergeometric_probabilities_b(sample_size, p_tree->number_of_leaves());
    this->compute_all_root_path_costs(std::back_inserter(root_path_costs));
    subtree_size_probabilities.assign(p_tree->number_of_leaves(), Protected_number_type(Number_type(0.0)));
    node_probabilities.assign(p_tree->number_of_nodes(), Protected_number_type(Number_type(0.0)));
    costs_kth_exponent.assign(p_tree->number_of_nodes(), Protected_number_type(Number_type(1.0)));

    std::set<int> unique_subtree_sizes;

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      unique_subtree_sizes.insert(p_tree->node(i).all_subtree_leaves);

    typename std::set<int>::iterator sit;

    for( sit=unique_subtree_sizes.begin(); sit!=unique_subtree_sizes.end(); sit++)
    {
      subtree_size_probabilities[(*sit)-1] = this->compute_node_probability(*sit,sample_size);

      ////////////////////////////////////////////////////////////////////////////
      //if(this->compute_node_probability(*sit,sample_size) > Number_type(1.0)) //
      //  this->compute_node_probability(*sit,sample_size,false);               // 
      ////////////////////////////////////////////////////////////////////////////
    }

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      node_probabilities[i] = subtree_size_probabilities[p_tree->node(i).all_subtree_leaves-1];

    std::vector<Protected_number_type> final_probabilities = node_probabilities;

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      for( int j=0; j<p_tree->node(i).number_of_children(); j++ )
        final_probabilities[i] = final_probabilities[i] - node_probabilities[p_tree->node(i).children[j]];

    Protected_number_type mean;

    bool is_ultrametric = p_tree->is_ultrametric();

    for(int i=0; i<k; i++)
    {
      Protected_number_type ith_raw_moment(0.0);

      if(i==1)
        for(int j=0; j < p_tree->number_of_nodes(); j++)
          root_path_costs[j] = root_path_costs[j]-mean;

      for(int j=0; j < p_tree->number_of_nodes(); j++)
      {
        costs_kth_exponent[j] = root_path_costs[j].pos_power(i+1);
        ith_raw_moment += costs_kth_exponent[j]*final_probabilities[j];
      }

      // The next check fixes a case where a floating point error occurs.

      if(sample_size==1 && i>0 && is_ultrametric == true)
       *ot++ = Number_type(0.0);
      else
       *ot++ = ith_raw_moment.to_number_type();

      if(i==0)
        mean = ith_raw_moment;

    } // for(int i=0; i<k; i++)

  } //compute_first_k_centralised_moments(int k, int sample_size, OutputIterator ot)

  template< class KernelType >
  template<class OutputIterator>
  void PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_first_k_raw_moments(int k, int sample_size, OutputIterator ot)
  {

    std::vector<Protected_number_type> temp;

    compute_first_k_raw_moments_protected(k, sample_size, std::back_inserter(temp));
    
    for(int i=0; i<temp.size(); i++)
      *ot++ = temp[i].to_number_type();

  } // compute_first_k_raw_moments(...)

#endif //CORE_ANCESTOR_COST_IMPL_H
