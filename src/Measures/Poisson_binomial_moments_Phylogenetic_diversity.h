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

#ifndef POISSON_BINOMIAL_MOMENTS_PHYLOGENETIC_DIVERSITY
#define POISSON_BINOMIAL_MOMENTS_PHYLOGENETIC_DIVERSITY

#include <thread>

namespace PhylogeneticMeasures
{

  template <class KernelType>
  class Poisson_binomial_moments_Phylogenetic_diversity
  {
   public:
  
    typedef KernelType                                   Kernel;
    typedef typename Kernel::Number_type                Number_type;
    typedef typename Kernel::Protected_number_type      Protected_number_type;
    typedef typename Kernel::Numeric_traits             Numeric_traits;
    typedef typename Kernel::Edge_relation_type         Edge_relation_type;
    typedef typename Numeric_traits::Square_root        Square_root;
    typedef typename Kernel::Polynomial                 Polynomial;
    typedef typename Kernel::Polynomial_multiplication  Polynomial_multiplication;
    typedef typename Kernel::Unimodal_tree              Tree_type;    
    typedef typename Tree_type::Node_type               Node_type; 
    typedef typename Tree_type::Leaves_iterator         Leaves_iterator;     
  
    typedef typename Kernel::Exception_type             Exception_type;
    typedef typename Kernel::Exception_functor          Exception_functor;

  private: 

    struct Node_polynomials
    {
      Polynomial basic, aggregated, aggregated_var;

      std::vector<Protected_number_type> exp_edge_contr_sums, var_edge_contr_sums;

      void clear()
      {
        basic.clear();
        aggregated.clear();
        aggregated_var.clear();
 
        exp_edge_contr_sums.clear();
        var_edge_contr_sums.clear();
      } 

    }; // struct Node_polynomials

   private:

    void _construct_node_levels_recursive(Tree_type &tree, int index, int level)
    {
      if( level < _node_levels.size())
        _node_levels[level].push_back(index);
      else
      {
        _node_levels.push_back(std::vector<int>());
        _node_levels[level].push_back(index);
      }

      Node_type v = tree.node(index);

      for(int i=0; i<v.number_of_children(); i++)
        _construct_node_levels_recursive(tree, v.children[i], level+1);
    }

    void _compute_polynomials_level_based(Tree_type &tree, int sample_size, 
                                          Polynomial &basic, Polynomial &aggregated, Polynomial &aggregated_var,
                                          std::vector<Protected_number_type> &exp_edge_contr_sums, 
                                          std::vector<Protected_number_type> &var_edge_contr_sums)
    {
      construct_node_levels(tree);

      _node_polynomials.assign(tree.number_of_nodes(), Node_polynomials());

      for(int i=_node_levels.size()-1; i>=0; i--)
      {
        for(int j=0; j<_node_levels[i].size(); j++)
        { 
          int index = _node_levels[i][j];

          _compute_polynomials_recursive(tree, index, sample_size, 
          _node_polynomials[index].basic, _node_polynomials[index].aggregated, 
          _node_polynomials[index].aggregated_var, _node_polynomials[index].exp_edge_contr_sums, 
          _node_polynomials[index].var_edge_contr_sums, true);

          Node_type v = tree.node(index);

          for(int k=0; k<v.number_of_children(); k++)
            _node_polynomials[v.children[k]].clear();
        }  

      } // for(int i=0; i<_node_levels.size(); i++)

      basic = _node_polynomials[tree.root_index()].basic;
      aggregated = _node_polynomials[tree.root_index()].aggregated;
      aggregated_var = _node_polynomials[tree.root_index()].aggregated_var;

      exp_edge_contr_sums = _node_polynomials[tree.root_index()].exp_edge_contr_sums;
      var_edge_contr_sums = _node_polynomials[tree.root_index()].var_edge_contr_sums;

    } // _compute_polynomials_level_based(...)
  
   public:
  
    Poisson_binomial_moments_Phylogenetic_diversity(): 
    _parallel_multiplication(true), _threads_mult(std::max(int(std::thread::hardware_concurrency()),1)){}

    void set_execution_parameters( bool par_mult, int threads_mult)
    {
      _parallel_multiplication = par_mult;
      _threads_mult = threads_mult;
    }

    void construct_node_levels(Tree_type &tree)
    { _construct_node_levels_recursive(tree,tree.root_index(),0);}
  
    template< class OutputIterator>
    void compute_expectations_and_variances
    ( Tree_type &in_tree, int sample_size, 
      OutputIterator ot_exp, OutputIterator ot_var,
      bool level_based = true)  
    {  
      Polynomial basic, aggregated, aggregated_var;
  
      Tree_type tree = in_tree;

      tree.convert_to_binary_tree();
      tree.assign_all_subtree_leaves();

      std::vector<Protected_number_type> exp_edge_contr_sums, var_edge_contr_sums; 

      if(level_based == false)
        _compute_polynomials_recursive( tree, tree.root_index(), 
                                        sample_size, basic, aggregated, 
                                        aggregated_var,exp_edge_contr_sums, 
                                        var_edge_contr_sums,false); 
      else  
        _compute_polynomials_level_based( tree, sample_size, basic, aggregated, 
                                          aggregated_var,exp_edge_contr_sums, 
                                          var_edge_contr_sums);  

      for(int i=0; i<basic.size(); i++)
      {
        Protected_number_type pn = (aggregated[i]-exp_edge_contr_sums[i])/basic[i],
                              var_minus = var_edge_contr_sums[i]/basic[i];

        *ot_exp++ = pn.to_number_type();

        pn = (aggregated_var[i]/basic[i])-var_minus-(pn*pn); 

        if(pn.to_number_type() < Number_type(0.0)) 
          *ot_var++ = Number_type(0.0);
        else
          *ot_var++ = pn.to_number_type();

      }

      if(sample_size > basic.size()-1)
        for(int i=basic.size(); i<=sample_size; i++)
        {
          *ot_exp++ = Number_type(0.0);
          *ot_var++ = Number_type(0.0);
        }

    } // compute_expectations_and_variances(...) 
  
   private:
  

    void _compute_polynomials_recursive
    ( Tree_type &tree, int index, int sample_size, 
      Polynomial &basic, Polynomial &aggregated, Polynomial &aggregated_var, 
      std::vector<Protected_number_type> &exp_edge_contr_sums, 
      std::vector<Protected_number_type> &var_edge_contr_sums,
      bool level_based = false)
    {  
      Node_type v = tree.node(index);
      Number_type distance(v.distance);
  
      if(v.number_of_children() == 0)
      {
        basic.push_back(Number_type(Number_type(1.0) - Number_type(tree.node_probability(index)))); 

        if(tree.node_probability(index) != Number_type(0.0))
          basic.push_back(tree.node_probability(index)); 

        Protected_number_type contr_0(Number_type(0.0)),
                              contr_1(distance);

        aggregated.push_back(basic[0]*contr_0);

        if(tree.node_probability(index) != Number_type(0.0)) 
          aggregated.push_back(basic[1]*contr_1);

        aggregated_var.push_back(basic[0]*contr_0*contr_0);

        if(tree.node_probability(index) != Number_type(0.0)) 
          aggregated_var.push_back(basic[1]*contr_1*contr_1);

        ////////////////////////////////////////////
        // Contributions that will be subtracted  //     
        ////////////////////////////////////////////

        exp_edge_contr_sums.push_back(basic[0]*contr_0);

        if(tree.node_probability(index) != Number_type(0.0)) 
          exp_edge_contr_sums.push_back(basic[1]*contr_1);
        else
          exp_edge_contr_sums.push_back(Protected_number_type(0.0));

        var_edge_contr_sums.push_back(basic[0]*contr_0*contr_0);

        if(tree.node_probability(index) != Number_type(0.0)) 
          var_edge_contr_sums.push_back(basic[1]*contr_1*contr_1);
        else
          var_edge_contr_sums.push_back(Protected_number_type(0.0));

      } // if(v.number_of_children() == 0)
      else
      {
        Polynomial_multiplication mult;

        mult.set_cutoff_value(sample_size);
        mult.set_is_parallel(_parallel_multiplication);
        mult.set_max_threads(_threads_mult);
              
        std::vector<Polynomial> ch_basic, ch_aggregated, ch_aggregated_var;
        std::vector< std::vector<Protected_number_type> > ch_exp_edge_contr_sums, 
                                                          ch_var_edge_contr_sums;  

        ch_basic.assign(v.number_of_children(), Polynomial());
        ch_aggregated.assign(v.number_of_children(), Polynomial());
        ch_aggregated_var.assign(v.number_of_children(), Polynomial());
        ch_exp_edge_contr_sums.assign(v.number_of_children(), Polynomial());
        ch_var_edge_contr_sums.assign(v.number_of_children(), Polynomial());

        if(level_based == true)
        {
          for(int i=0; i<v.number_of_children(); i++)
          {
            ch_basic[i] = _node_polynomials[v.children[i]].basic;
            ch_aggregated[i] = _node_polynomials[v.children[i]].aggregated;
            ch_aggregated_var[i] = _node_polynomials[v.children[i]].aggregated_var;

            ch_exp_edge_contr_sums[i] = _node_polynomials[v.children[i]].exp_edge_contr_sums;
            ch_var_edge_contr_sums[i] = _node_polynomials[v.children[i]].var_edge_contr_sums;
          }  

        } // if(level_based == true)
        else    
          for(int i=0; i<v.number_of_children(); i++)
            _compute_polynomials_recursive(tree,v.children[i],sample_size, 
                                           ch_basic[i], ch_aggregated[i], ch_aggregated_var[i], 
                                           ch_exp_edge_contr_sums[i], ch_var_edge_contr_sums[i],false);

        exp_edge_contr_sums.assign(v.all_subtree_leaves+1, Protected_number_type(0.0));    
        var_edge_contr_sums.assign(v.all_subtree_leaves+1, Protected_number_type(0.0));    
  
        // TODO: Check if it pays off in practice to multiply these 
        // polynomials in order of increasing size.
  
        for(int i=0; i<v.number_of_children(); i++)
        {
          if(i==0)
            basic = ch_basic[i];
          else
          {
            Polynomial tmp;
  
            mult(basic,ch_basic[i],tmp);
            basic = tmp;
          }
  
          Polynomial aggr, tmp_aggr, tmp_aggr_var;
  
          //////////////////////////////////////////////////////
          // Compute the product of all basic polynomials     //  
          // of the children, except for the one of child i.  //
          //////////////////////////////////////////////////////
  
          for(int j=0; j<v.number_of_children(); j++)
            if(j != i)
            {
              if(aggr.size() == 0)
                aggr = ch_basic[j]; 
              else
              {
                Polynomial tmp;
  
                mult(aggr,ch_basic[j],tmp);
                aggr = tmp;
              }
  
            } // if(j != i)
    
          //////////////////////////////////////////////////////////////////////
          // Multiply the resulting product of all basic polynomials          //
          // but the i-th with the aggregated polynomial of the i-th child.   //
          // Add the result to the aggregated polynomial of the current node. //
          // Do the same process also for the aggregated polynomial which is  //
          // used for the variance computations.                              //   
          //////////////////////////////////////////////////////////////////////
  
          if(aggregated.size() == 0)
            mult(ch_aggregated[i], aggr, aggregated);
          else
          {
            mult(ch_aggregated[i], aggr, tmp_aggr);
  
            for(int k=0; k<aggregated.size(); k++)
              if(k<tmp_aggr.size())
                aggregated[k] = aggregated[k] + tmp_aggr[k];
  
            for(int k=aggregated.size(); k<tmp_aggr.size(); k++)
              aggregated.push_back(tmp_aggr[k]);
          }

          if(aggregated_var.size() == 0)
            mult(ch_aggregated_var[i], aggr, aggregated_var);
          else
          { 
            mult(ch_aggregated_var[i], aggr, tmp_aggr_var);
   
            for(int k=0; k<aggregated_var.size(); k++)
              if(k<tmp_aggr_var.size())
                aggregated_var[k] = aggregated_var[k] + tmp_aggr_var[k];
  
            for(int k=aggregated_var.size(); k<tmp_aggr_var.size(); k++)
              aggregated_var.push_back(tmp_aggr_var[k]);
          }

          /////////////////////////////////////////////
          // Contributions that will be subtracted   //
          /////////////////////////////////////////////

          for(int tt=0; tt< ch_exp_edge_contr_sums[i].size(); tt++)
          {
            ch_exp_edge_contr_sums[i][tt] = ch_exp_edge_contr_sums[i][tt] * aggr[0];
            exp_edge_contr_sums[tt] +=  ch_exp_edge_contr_sums[i][tt];
          }

          for(int tt=0; tt< ch_var_edge_contr_sums[i].size(); tt++)
          {
            ch_var_edge_contr_sums[i][tt] = ch_var_edge_contr_sums[i][tt] * aggr[0];
            var_edge_contr_sums[tt] +=  ch_var_edge_contr_sums[i][tt];
          }

        ///////////////////////////////////////////////////////
        // For the variance polynomial, compute              //
        // the product of the aggregated polynomials of each //
        // pair of children separately, and add each product //
        // to aggregated_var.                                // 
        ///////////////////////////////////////////////////////


          for(int j=0; j<v.number_of_children(); j++)
            if(j < i) // This time we pick each possible pair once
            {
              Polynomial tmp;
  
              mult(ch_aggregated[i],ch_aggregated[j],tmp);

              for(int h=0; h<v.number_of_children(); h++)
                if( h!=i && h!=j )
                {
                  Polynomial sec_tmp;

                  mult(tmp,ch_basic[h],sec_tmp);
                  tmp = sec_tmp;
                } 

              for(int kk=0; kk<tmp.size(); kk++)
              {
                Protected_number_type pn = Protected_number_type(Number_type(2.0))*tmp[kk];

                if(kk<aggregated_var.size())
                  aggregated_var[kk] = aggregated_var[kk] + pn;
                else
                  aggregated_var.push_back(pn);
              }
  
            } // if(j < i)
  
        } // for(int i=0; i<v.number_of_children(); i++)


        /////////////////////////////////////////////
        // Contributions that will be subtracted   //
        /////////////////////////////////////////////

        if( tree.root_index() != index )
        {
          for(int tt=0; tt<basic.size(); tt++)
          {
            Protected_number_type contr;

            if(tt==0 || tt==sample_size) 
              contr = Protected_number_type(Number_type(0.0));
            else 
              contr = Protected_number_type(distance);

            exp_edge_contr_sums[tt] = exp_edge_contr_sums[tt] + (contr*basic[tt]);
          }

        } // if( tree.root_index() != index )

        //////////////////////////////////////////////
        // For the variance polynomial, compute     //
        // the two-times dot product between the    //
        // contribution of the current edge and     //
        // the contribution of all offspring edges. //
        // NOTE: to do this we need that polynomial //
        // 'aggregated' does not yet include the    //
        // contribution of the current edge.        // 
        //                                          //
        // We also add the squared contribution     //
        // of the current edge explicitly.          //
        //////////////////////////////////////////////


        if( tree.root_index() != index )
          for(int i=0; i<aggregated.size(); i++)
          {
            Protected_number_type contr; 

            if(i==0 || i==sample_size)
              contr = Protected_number_type(Number_type(0.0));
            else
              contr = Protected_number_type(distance);

            Protected_number_type pn = aggregated[i]*Protected_number_type(Number_type(2.0))*contr;

            pn += basic[i]*Protected_number_type(contr*contr);
            var_edge_contr_sums[i] = var_edge_contr_sums[i] + pn;

            if(i<aggregated_var.size())
              aggregated_var[i] = aggregated_var[i] + pn;
            else
              aggregated_var.push_back(pn); 
          }

        /////////////////////////////////////////  
        // Add the weighted polynomial of the  //
        // current node to the aggregated one. //
        /////////////////////////////////////////

        if( tree.root_index() != index )
          for(int i=0; i<basic.size(); i++)
          {
            Protected_number_type contr;

            if(i==0 || i== sample_size)
              contr = Protected_number_type(Number_type(0.0));
            else
              contr = Protected_number_type(distance);

            Protected_number_type pn = basic[i]*contr;

            if(i<aggregated.size())
              aggregated[i] = aggregated[i] + pn;
            else
              aggregated.push_back(pn);
          }


        ///////////////////////////////////////////////////
        // Kick out the coefficients that are not needed //
        ///////////////////////////////////////////////////

        while(aggregated.size()>sample_size+1)
          aggregated.pop_back();

        while(basic.size()>sample_size+1)
          basic.pop_back();

        while(aggregated_var.size()>sample_size+1)
          aggregated_var.pop_back();

      } // else of if(v.number_of_children == 0)

    } // _compute_polynomials_recursive(...)

   private:

    std::vector<std::vector<int> > _node_levels;
    std::vector< Node_polynomials > _node_polynomials;

    bool _parallel_multiplication; 
    int _threads_mult;
   
  }; // class Poisson_binomial_moments_Phylogenetic_diversity

} // namespace PhylogeneticMeasures

#endif // POISSON_BINOMIAL_MOMENTS_PHYLOGENETIC_DIVERSITY
