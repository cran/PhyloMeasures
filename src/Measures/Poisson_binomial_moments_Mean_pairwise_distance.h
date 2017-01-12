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

#ifndef POISSON_BINOMIAL_MOMENTS_MEAN_PAIRWISE_DISTANCE
#define POISSON_BINOMIAL_MOMENTS_MEAN_PAIRWISE_DISTANCE

#include <thread>

namespace PhylogeneticMeasures
{

  template <class KernelType>
  class Poisson_binomial_moments_Mean_pairwise_distance
  {
   public:
  
    typedef KernelType                                               Kernel;
    typedef Poisson_binomial_moments_Mean_pairwise_distance<Kernel>  Self;
    typedef typename Kernel::Number_type                            Number_type;
    typedef typename Kernel::Protected_number_type                  Protected_number_type;
    typedef typename Kernel::Numeric_traits                         Numeric_traits;
    typedef typename Kernel::Edge_relation_type                     Edge_relation_type;
    typedef typename Numeric_traits::Square_root                    Square_root;
    typedef typename Kernel::Polynomial                             Polynomial;
    typedef typename Kernel::Polynomial_multiplication              Polynomial_multiplication;
    typedef typename Kernel::Unimodal_tree                          Tree_type;    
    typedef typename Tree_type::Node_type                           Node_type;    
    typedef typename Tree_type::Leaves_iterator                     Leaves_iterator;    
    typedef typename Kernel::Exception_type                         Exception_type;
    typedef typename Kernel::Exception_functor                      Exception_functor;
  

   private: 

    struct Node_polynomials
    {
      Polynomial basic, poly_K, poly_KS, poly_KK, poly_KKS_KSK, poly_KSKS;

      std::vector<Protected_number_type> contr_sums_K, contr_sums_KS, 
                                         contr_sums_KK, contr_sums_KKS_KSK, contr_sums_KSKS;   

      void clear()
      {
        basic.clear();
        poly_K.clear();
        poly_KS.clear();
        poly_KK.clear();
        poly_KKS_KSK.clear();
        poly_KSKS.clear();
 
        contr_sums_K.clear();
        contr_sums_KS.clear();
        contr_sums_KK.clear();
        contr_sums_KKS_KSK.clear();
        contr_sums_KSKS.clear();
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
                                          Polynomial &basic, Polynomial &poly_K, Polynomial &poly_KS,
                                          Polynomial &poly_KK, Polynomial &poly_KKS_KSK, Polynomial &poly_KSKS,
                                          std::vector<Protected_number_type> &contr_sums_K, 
                                          std::vector<Protected_number_type> &contr_sums_KS, 
                                          std::vector<Protected_number_type> &contr_sums_KK, 
                                          std::vector<Protected_number_type> &contr_sums_KKS_KSK, 
                                          std::vector<Protected_number_type> &contr_sums_KSKS)
    {
      construct_node_levels(tree);

      _node_polynomials.assign(tree.number_of_nodes(), Node_polynomials());

      for(int i=_node_levels.size()-1; i>=0; i--)
      {
        for(int j=0; j<_node_levels[i].size(); j++)
        { 
          int index = _node_levels[i][j];

          _compute_polynomials_recursive(tree, index, sample_size, 
          _node_polynomials[index].basic, _node_polynomials[index].poly_K, 
          _node_polynomials[index].poly_KS, _node_polynomials[index].poly_KK, 
          _node_polynomials[index].poly_KKS_KSK, _node_polynomials[index].poly_KSKS,
          _node_polynomials[index].contr_sums_K, _node_polynomials[index].contr_sums_KS,
          _node_polynomials[index].contr_sums_KK, _node_polynomials[index].contr_sums_KKS_KSK,
          _node_polynomials[index].contr_sums_KSKS, true, &_node_polynomials);

          Node_type v = tree.node(index);

          for(int k=0; k<v.number_of_children(); k++)
            _node_polynomials[v.children[k]].clear();
        }  

      } // for(int i=0; i<_node_levels.size(); i++)


      basic = _node_polynomials[tree.root_index()].basic;
      poly_K = _node_polynomials[tree.root_index()].poly_K;
      poly_KS = _node_polynomials[tree.root_index()].poly_KS;
      poly_KK = _node_polynomials[tree.root_index()].poly_KK;
      poly_KKS_KSK = _node_polynomials[tree.root_index()].poly_KKS_KSK;
      poly_KSKS = _node_polynomials[tree.root_index()].poly_KSKS;

      contr_sums_K = _node_polynomials[tree.root_index()].contr_sums_K;
      contr_sums_KS = _node_polynomials[tree.root_index()].contr_sums_KS;
      contr_sums_KK = _node_polynomials[tree.root_index()].contr_sums_KK;
      contr_sums_KKS_KSK = _node_polynomials[tree.root_index()].contr_sums_KKS_KSK;
      contr_sums_KSKS = _node_polynomials[tree.root_index()].contr_sums_KSKS;

    } // _compute_polynomials_level_based(...)

   public:
  
    Poisson_binomial_moments_Mean_pairwise_distance(): 
    _parallel_multiplication(true), _threads_mult(std::max(int(std::thread::hardware_concurrency()),1)){}

    void set_execution_parameters( bool par_mult, int threads_mult)
    {
      _parallel_multiplication = par_mult;
      _threads_mult = threads_mult;
    }

    void construct_node_levels(Tree_type &tree)
    { _construct_node_levels_recursive(tree,tree.root_index(),0);}
  
    template<class OutputIterator>
    void compute_expectations_and_variances
    ( Tree_type &in_tree, int sample_size, 
      OutputIterator ot_exp, OutputIterator ot_var, 
      bool level_based = true )  
    {
      if(sample_size >=0 )
      {
        *ot_exp++ = Number_type(0.0);
        *ot_var++ = Number_type(0.0);
      }

      if(sample_size >=1 )
      {
        *ot_exp++ = Number_type(0.0);
        *ot_var++ = Number_type(0.0);
      }

      if(sample_size < 2)
        return;

      Tree_type tree = in_tree;

      tree.convert_to_binary_tree();
      tree.assign_all_subtree_leaves();

      Polynomial basic, poly_K, poly_KS, poly_KK, poly_KKS_KSK, poly_KSKS;
  
      std::vector<Protected_number_type> contr_sums_K, contr_sums_KS, 
                                         contr_sums_KK, contr_sums_KKS_KSK, contr_sums_KSKS; 

      if(level_based == false)
      {
        _compute_polynomials_recursive( tree, tree.root_index(), 
                                        sample_size, basic, poly_K, poly_KS,
                                        poly_KK, poly_KKS_KSK, poly_KSKS,
                                        contr_sums_K, contr_sums_KS, 
                                        contr_sums_KK, contr_sums_KKS_KSK, 
                                        contr_sums_KSKS, false); 
      }
      else
      {
        _compute_polynomials_level_based( tree, sample_size, basic, poly_K, poly_KS,
                                          poly_KK, poly_KKS_KSK, poly_KSKS,
                                          contr_sums_K, contr_sums_KS, 
                                          contr_sums_KK, contr_sums_KKS_KSK, 
                                          contr_sums_KSKS);
      }

      for(int i=2; i<basic.size(); i++)
      {
        Number_type ci(i);
        Number_type factor = Number_type(2.0)/(ci*(ci-Number_type(1.0)));
        Protected_number_type ss(ci), exp, var;
      
        exp = ((ss*(poly_K[i])) - (poly_KS[i]))/basic[i];     

        *ot_exp++ = factor*exp.to_number_type();

        var = (ss*ss*(poly_KK[i]-contr_sums_KK[i])) - 
              (ss*(poly_KKS_KSK[i]-contr_sums_KKS_KSK[i])) + 
              poly_KSKS[i]-contr_sums_KSKS[i];

        var = Protected_number_type(factor*factor)*((var/basic[i])-(exp*exp)); 

        if(var.to_number_type() < Number_type(0.0)) 
          *ot_var++ = Number_type(0.0);
        else
          *ot_var++ = var.to_number_type();
      }

      if(sample_size > basic.size()-1)
        for(int i=basic.size(); i<=sample_size; i++)
        {
          *ot_exp++ = Number_type(0.0);
          *ot_var++ = Number_type(0.0);
        }

    } // compute_expectations_and_variances(...) 
  
   protected:
    
    void _compute_polynomials_recursive
    ( Tree_type &tree, int index, int sample_size, 
      Polynomial &basic, Polynomial &poly_K, Polynomial &poly_KS,
      Polynomial &poly_KK, Polynomial &poly_KKS_KSK, Polynomial &poly_KSKS,
      std::vector<Protected_number_type> &contr_sums_K, 
      std::vector<Protected_number_type> &contr_sums_KS,
      std::vector<Protected_number_type> &contr_sums_KK,
      std::vector<Protected_number_type> &contr_sums_KKS_KSK,
      std::vector<Protected_number_type> &contr_sums_KSKS, 
      bool level_based = false, 
      std::vector< Node_polynomials > *node_polynomials = NULL)
    {  
      Node_type v = tree.node(index);
      Number_type distance(v.distance);
  
      if(v.number_of_children() == 0)
      {
        basic.push_back(Number_type(Number_type(1.0) - Number_type(tree.node_probability(index)))); 

        if(tree.node_probability(index) != Number_type(0.0))
          basic.push_back(tree.node_probability(index)); 

        Protected_number_type contr_0(Number_type(0.0)),
                              contr_1(Number_type(1.0)*distance);

        poly_K.push_back(contr_0);

        if(tree.node_probability(index) != Number_type(0.0)) 
          poly_K.push_back(basic[1]*contr_1);

        poly_KS.push_back(contr_0);

        if(tree.node_probability(index) != Number_type(0.0)) 
          poly_KS.push_back(basic[1]*contr_1);

        poly_KK.push_back(contr_0);

        if(tree.node_probability(index) != Number_type(0.0)) 
          poly_KK.push_back(basic[1]*contr_1*contr_1);

        poly_KKS_KSK.push_back(Protected_number_type(Number_type(0.0)));

        if(tree.node_probability(index) != Number_type(0.0)) 
          poly_KKS_KSK.push_back(Protected_number_type(Number_type(2.0))*basic[1]*contr_1*contr_1);

        poly_KSKS.push_back(Protected_number_type(Number_type(0.0)));

        if(tree.node_probability(index) != Number_type(0.0)) 
          poly_KSKS.push_back(basic[1]*contr_1*contr_1);

        ////////////////////////////////////////////
        // Contributions that will be subtracted  //     
        ////////////////////////////////////////////

        contr_sums_K.push_back(contr_0);

        if(tree.node_probability(index) != Number_type(0.0)) 
          contr_sums_K.push_back(basic[1]*contr_1);
        else
          contr_sums_K.push_back(Protected_number_type(0.0));


        contr_sums_KS.push_back(contr_0);

        if(tree.node_probability(index) != Number_type(0.0)) 
          contr_sums_KS.push_back(basic[1]*contr_1);
        else
          contr_sums_KS.push_back(Protected_number_type(0.0));

        contr_sums_KK.push_back(contr_0);

        if(tree.node_probability(index) != Number_type(0.0)) 
          contr_sums_KK.push_back(basic[1]*contr_1*contr_1);
        else
          contr_sums_KK.push_back(Protected_number_type(0.0));


        contr_sums_KKS_KSK.push_back(Protected_number_type(Number_type(0.0)));

        if(tree.node_probability(index) != Number_type(0.0)) 
          contr_sums_KKS_KSK.push_back(Protected_number_type(Number_type(2.0))*basic[1]*contr_1*contr_1);
        else
          contr_sums_KKS_KSK.push_back(Protected_number_type(Number_type(0.0)));

        contr_sums_KSKS.push_back(Protected_number_type(Number_type(0.0)));

        if(tree.node_probability(index) != Number_type(0.0)) 
          contr_sums_KSKS.push_back(basic[1]*contr_1*contr_1);
        else
          contr_sums_KSKS.push_back(Protected_number_type(Number_type(0.0)));

      } // if(v.number_of_children() == 0)
      else
      {              
        Polynomial_multiplication mult;

        mult.set_cutoff_value(sample_size);
        mult.set_is_parallel(_parallel_multiplication);
        mult.set_max_threads(_threads_mult); 

        std::vector<Polynomial> ch_basic, ch_poly_K, ch_poly_KK, 
                                ch_poly_KS, ch_poly_KKS_KSK, ch_poly_KSKS;

        std::vector< std::vector<Protected_number_type> > ch_contr_sums_K, ch_contr_sums_KK,
                                                          ch_contr_sums_KS, ch_contr_sums_KKS_KSK,
                                                          ch_contr_sums_KSKS;
  
        ch_basic.assign(v.number_of_children(), Polynomial());
        ch_poly_K.assign(v.number_of_children(), Polynomial());
        ch_poly_KK.assign(v.number_of_children(), Polynomial());
        ch_poly_KS.assign(v.number_of_children(), Polynomial());
        ch_poly_KKS_KSK.assign(v.number_of_children(), Polynomial());
        ch_poly_KSKS.assign(v.number_of_children(), Polynomial());

        ch_contr_sums_K.assign(v.number_of_children(), Polynomial());
        ch_contr_sums_KK.assign(v.number_of_children(), Polynomial());
        ch_contr_sums_KS.assign(v.number_of_children(), Polynomial());
        ch_contr_sums_KKS_KSK.assign(v.number_of_children(), Polynomial());   
        ch_contr_sums_KSKS.assign(v.number_of_children(), Polynomial());


        if(level_based == true)
        {
          for(int i=0; i<v.number_of_children(); i++)
          {
            ch_basic[i] = (*node_polynomials)[v.children[i]].basic;
            ch_poly_K[i] = (*node_polynomials)[v.children[i]].poly_K;
            ch_poly_KK[i] = (*node_polynomials)[v.children[i]].poly_KK;
            ch_poly_KS[i] = (*node_polynomials)[v.children[i]].poly_KS;
            ch_poly_KKS_KSK[i] = (*node_polynomials)[v.children[i]].poly_KKS_KSK;
            ch_poly_KSKS[i] = (*node_polynomials)[v.children[i]].poly_KSKS;

            ch_contr_sums_K[i] = (*node_polynomials)[v.children[i]].contr_sums_K;
            ch_contr_sums_KK[i] = (*node_polynomials)[v.children[i]].contr_sums_KK;
            ch_contr_sums_KS[i] = (*node_polynomials)[v.children[i]].contr_sums_KS;
            ch_contr_sums_KKS_KSK[i] = (*node_polynomials)[v.children[i]].contr_sums_KKS_KSK;
            ch_contr_sums_KSKS[i] = (*node_polynomials)[v.children[i]].contr_sums_KSKS;
          }  

        } // if(level_based == true)
        else 
          for(int i=0; i<v.number_of_children(); i++)
            _compute_polynomials_recursive(tree,v.children[i],sample_size, 
                                           ch_basic[i], ch_poly_K[i], ch_poly_KS[i],
                                           ch_poly_KK[i], ch_poly_KKS_KSK[i], ch_poly_KSKS[i], 
                                           ch_contr_sums_K[i], ch_contr_sums_KS[i],
                                           ch_contr_sums_KK[i], ch_contr_sums_KKS_KSK[i],
                                           ch_contr_sums_KSKS[i],false);

        contr_sums_K.assign(v.all_subtree_leaves+1, Protected_number_type(0.0));    
        contr_sums_KS.assign(v.all_subtree_leaves+1, Protected_number_type(0.0));    
        contr_sums_KK.assign(v.all_subtree_leaves+1, Protected_number_type(0.0));    
        contr_sums_KKS_KSK.assign(v.all_subtree_leaves+1, Protected_number_type(0.0));    
        contr_sums_KSKS.assign(v.all_subtree_leaves+1, Protected_number_type(0.0));    
 
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
          // but the i-th with the poly_K polynomial of the i-th child.       //
          // Add the result to the poly_K polynomial of the current node.     //
          // Do the same process also for the poly_K polynomial which is      //
          // used for the variance computations.                              //   
          //////////////////////////////////////////////////////////////////////
  
          if(poly_K.size() == 0)
            mult(ch_poly_K[i], aggr, poly_K);
          else
          {
            mult(ch_poly_K[i], aggr, tmp_aggr);
  
            for(int k=0; k<poly_K.size(); k++)
              if(k<tmp_aggr.size())
                poly_K[k] = poly_K[k] + tmp_aggr[k];
  
            for(int k=poly_K.size(); k<tmp_aggr.size(); k++)
              poly_K.push_back(tmp_aggr[k]);
          }

          tmp_aggr.clear();

          if(poly_KS.size() == 0)
            mult(ch_poly_KS[i], aggr, poly_KS);
          else
          {
            mult(ch_poly_KS[i], aggr, tmp_aggr);
  
            for(int k=0; k<poly_KS.size(); k++)
              if(k<tmp_aggr.size())
                poly_KS[k] = poly_KS[k] + tmp_aggr[k];
  
            for(int k=poly_KS.size(); k<tmp_aggr.size(); k++)
              poly_KS.push_back(tmp_aggr[k]);
          }

          if(poly_KK.size() == 0)
            mult(ch_poly_KK[i], aggr, poly_KK);
          else
          {
            mult(ch_poly_KK[i], aggr, tmp_aggr_var);
   
            for(int k=0; k<poly_KK.size(); k++)
              if(k<tmp_aggr_var.size())
                poly_KK[k] = poly_KK[k] + tmp_aggr_var[k];
  
            for(int k=poly_KK.size(); k<tmp_aggr_var.size(); k++)
              poly_KK.push_back(tmp_aggr_var[k]);
          }

          tmp_aggr_var.clear();

          if(poly_KSKS.size() == 0)
            mult(ch_poly_KSKS[i], aggr, poly_KSKS);
          else
          {
            mult(ch_poly_KSKS[i], aggr, tmp_aggr_var);
   
            for(int k=0; k<poly_KSKS.size(); k++)
              if(k<tmp_aggr_var.size())
                poly_KSKS[k] = poly_KSKS[k] + tmp_aggr_var[k];
  
            for(int k=poly_KSKS.size(); k<tmp_aggr_var.size(); k++)
              poly_KSKS.push_back(tmp_aggr_var[k]);
          }

          tmp_aggr_var.clear();

          if(poly_KKS_KSK.size() == 0)
            mult(ch_poly_KKS_KSK[i], aggr, poly_KKS_KSK);
          else
          {
            mult(ch_poly_KKS_KSK[i], aggr, tmp_aggr_var);
   
            for(int k=0; k<poly_KKS_KSK.size(); k++)
              if(k<tmp_aggr_var.size())
                poly_KKS_KSK[k] = poly_KKS_KSK[k] + tmp_aggr_var[k];
  
            for(int k=poly_KKS_KSK.size(); k<tmp_aggr_var.size(); k++)
              poly_KKS_KSK.push_back(tmp_aggr_var[k]);
          }

          /////////////////////////////////////////////
          // Contributions that will be subtracted   //
          /////////////////////////////////////////////

          for(int tt=0; tt< ch_contr_sums_K[i].size(); tt++)
          {
            ch_contr_sums_K[i][tt] = ch_contr_sums_K[i][tt] * aggr[0];
            contr_sums_K[tt] +=  ch_contr_sums_K[i][tt];
          }

          for(int tt=0; tt< ch_contr_sums_KS[i].size(); tt++)
          {
            ch_contr_sums_KS[i][tt] = ch_contr_sums_KS[i][tt] * aggr[0];
            contr_sums_KS[tt] +=  ch_contr_sums_KS[i][tt];
          }

          for(int tt=0; tt< ch_contr_sums_KK[i].size(); tt++)
          {
            ch_contr_sums_KK[i][tt] = ch_contr_sums_KK[i][tt] * aggr[0];
            contr_sums_KK[tt] +=  ch_contr_sums_KK[i][tt];
          }

          for(int tt=0; tt< ch_contr_sums_KSKS[i].size(); tt++)
          {
            ch_contr_sums_KSKS[i][tt] = ch_contr_sums_KSKS[i][tt] * aggr[0];
            contr_sums_KSKS[tt] +=  ch_contr_sums_KSKS[i][tt];
          }

          for(int tt=0; tt< ch_contr_sums_KKS_KSK[i].size(); tt++)
          {
            ch_contr_sums_KKS_KSK[i][tt] = ch_contr_sums_KKS_KSK[i][tt] * aggr[0];
            contr_sums_KKS_KSK[tt] +=  ch_contr_sums_KKS_KSK[i][tt];
          }

          ///////////////////////////////////////////////////////
          // For the variance polynomial, compute              //
          // the product of the ch_poly_K polynomials of each  //
          // pair of children separately, and add each product //
          // to poly_KK.                                       // 
          ///////////////////////////////////////////////////////

          for(int j=0; j<v.number_of_children(); j++)
            if(j < i) // This time we pick each possible pair once
            {
              Polynomial tmp;
  
              mult(ch_poly_K[i],ch_poly_K[j],tmp);

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

                if(kk<poly_KK.size())
                  poly_KK[kk] = poly_KK[kk] + pn;
                else
                  poly_KK.push_back(pn); 
              }
  
            } // if(j < i)
 
          ///////////////////////////////////////////////////////
          // For the variance computations, compute            //
          // the product of the ch_poly_KS polynomials of each //
          // pair of children separately, and add each product //
          // to poly_KSKS.                                     // 
          ///////////////////////////////////////////////////////


          for(int j=0; j<v.number_of_children(); j++)
            if(j < i) // This time we pick each possible pair once
            {
              Polynomial tmp;
  
              mult(ch_poly_KS[i],ch_poly_KS[j],tmp);

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

                if(kk<poly_KSKS.size())
                  poly_KSKS[kk] = poly_KSKS[kk] + pn;
                else
                  poly_KSKS.push_back(pn); 
              }
  
            } // if(j < i)

          ///////////////////////////////////////////////////////
          // For the variance computations, compute            //
          // the product of the ch_poly_KS polynomials of each //
          // pair of children separately, and add each product //
          // to poly_KKS_KSK.                                  // 
          ///////////////////////////////////////////////////////

          for(int j=0; j<v.number_of_children(); j++)
            if(j < i) // This time we pick each possible pair once
            {
              Polynomial tmp1, tmp2;
  
              mult(ch_poly_KS[i],ch_poly_K[j],tmp1);
              mult(ch_poly_K[i],ch_poly_KS[j],tmp2);

              for(int h=0; h<v.number_of_children(); h++)
                if( h!=i && h!=j )
                {
                  Polynomial sec_tmp;

                  mult(tmp1,ch_basic[h],sec_tmp);
                  tmp1 = sec_tmp; 

                  sec_tmp.clear(); 

                  mult(tmp2,ch_basic[h],sec_tmp);
                  tmp2 = sec_tmp;                 
                } 

              for(int kk=0; kk<tmp1.size(); kk++)
              {
                Protected_number_type pn = Protected_number_type(Number_type(2.0))*tmp1[kk];

                if(kk<poly_KKS_KSK.size())
                  poly_KKS_KSK[kk] = poly_KKS_KSK[kk] + pn;
                else
                  poly_KKS_KSK.push_back(pn); 
              }

              for(int kk=0; kk<tmp2.size(); kk++)
              {
                Protected_number_type pn = Protected_number_type(Number_type(2.0))*tmp2[kk];

                if(kk<poly_KKS_KSK.size())
                  poly_KKS_KSK[kk] = poly_KKS_KSK[kk] + pn;
                else
                  poly_KKS_KSK.push_back(pn); 
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
            Number_type ntt(tt);
            Protected_number_type contr(ntt), dist(distance);
            contr_sums_K[tt] = contr_sums_K[tt] + (contr*dist*basic[tt]);
            contr_sums_KS[tt] = contr_sums_KS[tt] + (contr*contr*dist*basic[tt]);
          }

        } // if( tree.root_index() != index )

        ///////////////////////////////////////////////
        // For the variance computation, compute     //
        // the two-times dot product between the     //
        // contribution of the current edge and      //
        // the contribution of all offspring edges.  //
        // NOTE: to do this we need that polynomials //
        // 'poly_K' and 'poly_KS' do not yet include //
        // the contribution of the current edge.     // 
        //                                           //
        // We also add the squared contribution      //
        // of the current edge explicitly.           //
        ///////////////////////////////////////////////

        if( tree.root_index() != index )
          for(int i=0; i<poly_K.size(); i++)
          {
            Number_type contr(i);
            Protected_number_type pn = poly_K[i]*Protected_number_type(Number_type(2.0)*contr*distance);

            pn += basic[i]*Protected_number_type(contr*contr*distance*distance);

            contr_sums_KK[i] = contr_sums_KK[i] + pn;

            if(i<poly_KK.size())
              poly_KK[i] = poly_KK[i] + pn;
            else
              poly_KK.push_back(pn); 
          }

        if( tree.root_index() != index )
          for(int i=0; i<poly_KS.size(); i++)
          {
            Number_type contr(i);

            Protected_number_type tcc(Number_type(2.0)*contr*contr*distance),
                                  cccc(contr*contr*contr*contr*distance*distance);

            Protected_number_type pn = poly_KS[i]*tcc;

            pn = pn + (basic[i]*cccc);

            contr_sums_KSKS[i] = contr_sums_KSKS[i] + pn;

            if(i<poly_KSKS.size())
              poly_KSKS[i] = poly_KSKS[i] + pn;
            else
              poly_KSKS.push_back(pn);
          }


        if( tree.root_index() != index )
          for(int i=0; i<poly_KKS_KSK.size(); i++)
          {
            Number_type contr(i);
            Protected_number_type pn = poly_KS[i]*Protected_number_type(Number_type(2.0)*contr*distance);

            pn += poly_K[i]*Protected_number_type(Number_type(2.0)*contr*contr*distance);
            pn += basic[i]*Protected_number_type(Number_type(2.0)*contr*contr*contr*distance*distance);

            contr_sums_KKS_KSK[i] = contr_sums_KKS_KSK[i] + pn;

            if(i<poly_KKS_KSK.size())
              poly_KKS_KSK[i] = poly_KKS_KSK[i] + pn;
            else
              poly_KKS_KSK.push_back(pn); 
          }

        /////////////////////////////////////////  
        // Add the weighted polynomial of the  //
        // current node to poly_K.             //
        /////////////////////////////////////////
  
        if( tree.root_index() != index )
          for(int i=0; i<basic.size(); i++)
          {
            Number_type ni(i);
            Protected_number_type contr(ni);
            Protected_number_type pn = basic[i]*contr*distance, 
                                  pns = basic[i]*contr*contr*distance; 

            if(i<poly_K.size())
              poly_K[i] = poly_K[i] + pn;
            else
              poly_K.push_back(pn);


            if(i<poly_KS.size())
              poly_KS[i] = poly_KS[i] + pns;
            else
              poly_KS.push_back(pns); 

          } // for(int i=0; i<basic.size(); i++)


        ///////////////////////////////////////////////////
        // Kick out the coefficients that are not needed //
        ///////////////////////////////////////////////////

        //while(basic.size()>sample_size+1)
        //  basic.pop_back();

        //while(poly_K.size()>sample_size+1)
        //  poly_K.pop_back();

        //while(poly_KS.size()>sample_size+1)
        //  poly_KS.pop_back();

        //if(poly_KK.size()>sample_size+1)
        //  poly_KK.pop_back();

        //if(poly_KSKS.size()>sample_size+1)
        //  poly_KSKS.pop_back();

        //if(poly_KKS_KSK.size()>sample_size+1)
        //  poly_KKS_KSK.pop_back();

      } // else of if(v.number_of_children == 0)

    } // _compute_polynomials_recursive(...)

   private:

    std::vector<std::vector<int> > _node_levels;
    std::vector< Node_polynomials > _node_polynomials;

    bool _parallel_multiplication; 
    int  _threads_mult;
   
  }; // class Poisson_binomial_moments_Mean_pairwise_distance

} // namespace PhylogeneticMeasures

#endif // POISSON_BINOMIAL_MOMENTS_MEAN_PAIRWISE_DISTANCE
