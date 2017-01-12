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

#ifndef SEQUENTIAL_SAMPLER_H
#define SEQUENTIAL_SAMPLER_H

#include<vector>
#include<random>
#include<chrono>

namespace PhylogeneticMeasures{

template <class KernelType>
class Sequential_sampler
{

 public:

  typedef KernelType                               Kernel;
  typedef typename Kernel::Number_type            Number_type;
  typedef typename Kernel::Exception_type         Exception_type;
  typedef typename Kernel::Exception_functor      Exception_functor;

  struct Sampler_data
  {
    std::vector<int> *vals; 
    std::vector<Number_type> *abundances;
  };

 private:

  struct Sum_node_type 
  {
    Number_type abundance; 
    int lc, rc, parent, val, leaves_in_subtree;

    Sum_node_type():abundance(-1.0), lc(-1),rc(-1),parent(-1),val(-1), leaves_in_subtree(0){}

  }; // Sum_node_type

 private:

  int _select_random_element_recursive(Number_type p, int index)
  {
    Sum_node_type v = _tree[index];
  


    if(v.abundance == Number_type(0.0) || _n == 0 || v.leaves_in_subtree == 0)
    {
      std::string exception_msg;
      exception_msg += " Cannot sample more elements from the species pool.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(v.lc == -1 )
    {
      _sampled_nodes.push_back(index);
      _sampled_abundances.push_back(v.abundance);          
      _update_path_to_root(index, v.abundance); 
      return v.val;
    }

    // The two first conditions are provided for avoiding 
    // errors caused by numerical errors.

    if(_tree[v.rc].leaves_in_subtree == 0)
      return _select_random_element_recursive(p, v.lc);
    else if(_tree[v.lc].leaves_in_subtree == 0)
      return _select_random_element_recursive(p, v.rc);
    else if(p <= _tree[v.lc].abundance)
      return _select_random_element_recursive(p, v.lc);
    else
      return _select_random_element_recursive(p-_tree[v.lc].abundance, v.rc);      

  } // _select_random_element_recursive(...)

  void _update_path_to_root(int index, Number_type p)
  {
    Sum_node_type v = _tree[index];

    _tree[index].leaves_in_subtree--;
 
    if(v.leaves_in_subtree == 1)
      _tree[index].abundance = Number_type(0.0);
    else
      _tree[index].abundance = v.abundance - p;

    if(v.parent !=-1)
      _update_path_to_root(v.parent,p);

  } // _update_path_to_root(int index, Number_type p)

  void _restore(int index, Number_type abundance)
  {
    Sum_node_type v = _tree[index];

    _tree[index].leaves_in_subtree++;
    _tree[index].abundance += abundance;

    if(v.parent !=-1)
      _restore(v.parent,abundance);
  }

  void _initialize(std::vector<int> &vals, std::vector<Number_type> &abundances)
  {
    _n=vals.size();

    for(int i=0; i<vals.size(); i++)
      _vals.push_back(vals[i]);

    for(int i=0; i<abundances.size(); i++)
      _abundances.push_back(abundances[i]);

    if(abundances.size() != vals.size())
    {
      std::string exception_msg;
      exception_msg += " The number of species does not match";
      exception_msg += " the number of the provided abundance values.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }    

    for(int i=0; i<_n; i++)
    {
      Sum_node_type v;
  
      if(abundances[i] <= Number_type(0.0))
      {
        std::string exception_msg;
        exception_msg += " Negative or zero abundance values are not allowed.\n";     
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }       

      v.abundance = abundances[i];
      v.val = vals[i];
      v.leaves_in_subtree = 1;

      _tree.push_back(v);
    }

    _construct_tree();

  } // _initialize(...)

  int _construct_tree()
  {
    std::vector<int> nodes;   
    int last_index=_n-1;

    for(int i=0; i<_n; i++)
      nodes.push_back(i);

    while(nodes.size() > 1)
    {
      std::vector<int> new_nodes;

      for(int i=0; i<nodes.size(); i++)
        if(i%2 == 1)
        {
          Sum_node_type v;

          last_index++; 

          v.lc = nodes[i-1];
          v.rc = nodes[i];
          v.abundance = _tree[nodes[i-1]].abundance + _tree[nodes[i]].abundance;
          v.leaves_in_subtree = _tree[nodes[i-1]].leaves_in_subtree + _tree[nodes[i]].leaves_in_subtree;

          _tree[nodes[i-1]].parent = last_index;
          _tree[nodes[i]].parent = last_index;
            
          _tree.push_back(v);

          new_nodes.push_back(last_index);
        }   

      if(nodes.size()%2 == 1)
        new_nodes.push_back(nodes.back());

      nodes = new_nodes;

    } // while(nodes.size() > 1)

    return last_index;
  
  } // _construct_tree()

 public:

  Sequential_sampler(Sampler_data data, unsigned seed):
  _generator(seed),_distribution(0.0,1.0)
  { _initialize(*(data.vals),*(data.abundances)); }

  Sequential_sampler(Sampler_data data):
  _generator(std::chrono::system_clock::now().time_since_epoch().count()),_distribution(0.0,1.0)
  { _initialize(*(data.vals),*(data.abundances)); }

  Sequential_sampler(std::vector<int> &vals, std::vector<Number_type> &abundances, unsigned seed):
  _generator(seed),_distribution(0.0,1.0)
  { _initialize(vals,abundances); }

  Sequential_sampler(std::vector<int> &vals, std::vector<Number_type> &abundances):
  _generator(std::chrono::system_clock::now().time_since_epoch().count()),_distribution(0.0,1.0)
  { _initialize(vals,abundances); }

  Sampler_data data()
  {
    Sampler_data dt;

    dt.abundances = &_abundances;
    dt.vals = &_vals;

    return dt;
  }

  int select_random_element()
  {
    Number_type p=Number_type(_distribution(_generator))*_tree.back().abundance;

    return _select_random_element_recursive(p,_tree.size()-1);
  }

  void restore()
  {  
    for(int i=0; i<_sampled_nodes.size(); i++)
      _restore(_sampled_nodes[i], _sampled_abundances[i]); 

    _sampled_nodes.clear();
    _sampled_abundances.clear();

  } // restore()

  int root_index()
  { return _tree.size()-1;}

  void operator()(int sample_size, std::vector<int> &sample)
  {
    if(sample_size > _n || sample_size < 0)
    {
      std::string exception_msg;
      exception_msg += " Requested sample size is out of range.\n";
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    for(int i=0; i<sample_size; i++)
      sample.push_back(select_random_element());

    restore();

  } // operator()(...)

 private:

  int _n;
  std::vector<bool> _vb;
  std::vector<int> _vals, _sampled_nodes;
  std::vector<Number_type> _abundances, _sampled_abundances;
  std::vector<Sum_node_type> _tree;
  std::default_random_engine _generator;
  std::uniform_real_distribution<double> _distribution;

}; // class sequential_sampler

} // PhylogeneticMeasures

#endif // SEQUENTIAL_SAMPLER_H
