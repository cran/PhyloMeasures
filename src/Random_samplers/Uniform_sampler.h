//////////////////////////////////////////////////////////////////////////////////
//    Copyright (C) 2015,  Constantinos Tsirogiannis.  Email: tsirogiannis.c@gmail.com
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

#ifndef UNIFORM_SAMPLER_H
#define UNIFORM_SAMPLER_H

#include<vector>
#include<algorithm>
#include<random>
#include<chrono>

namespace PhylogeneticMeasures{

template <class KernelType>
class Uniform_sampler
{

 public:

  typedef KernelType                          Kernel;
  typedef typename Kernel::Number_type        Number_type;
  typedef typename Kernel::Exception_type     Exception_type;
  typedef typename Kernel::Exception_functor  Exception_functor;

  struct Sampler_data
  {
    std::vector<int> *vals;
  };

 private:

  void _initialize(std::vector<int> &vals)
  {
    for(int i=0; i<vals.size(); i++)
      if(vals[i]<0)
      {
        std::string exception_msg;
        exception_msg += " Uniform sampler can handle only non-negative values.\n";     
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }

    _n=vals.size();

    for(int i=0; i<vals.size(); i++)
      _vals.push_back(vals[i]);

  } // _initialize(std::vector<int> &vals)

 public:

  Uniform_sampler(Sampler_data data, unsigned int seed):
  _generator(seed), _distribution(0,data.vals->size()-1)
  { _initialize(*(data.vals));}

  Uniform_sampler(Sampler_data data):
  _generator(std::chrono::system_clock::now().time_since_epoch().count()), _distribution(0,data.vals->size()-1)
  { _initialize(*(data.vals));}

  Uniform_sampler(std::vector<int> &vals, unsigned int seed):
  _generator(seed), _distribution(0,vals.size()-1)
  { _initialize(vals);}

  Uniform_sampler(std::vector<int> &vals):
  _generator(std::chrono::system_clock::now().time_since_epoch().count()), _distribution(0,vals.size()-1)
  { _initialize(vals);}

  Sampler_data data()
  {
    Sampler_data dt;

    dt.vals = &_vals;

    return dt;
  }

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

    if(sample_size<=_n/2)
      select_random_sample(sample_size, sample);
    else
    {
      select_counter_sample(_n-sample_size, sample);
      std::shuffle(sample.begin(), sample.end(), _generator);
 
    } // else of if(sample_size<=n/2)

  } // operator()(...)

  void select_random_sample(int sample_size, std::vector<int> &sample)
  {
    int count=0;
    std::vector<int> indices;   

    while(count < sample_size)
    {
      int index = _distribution(_generator);

      if(_vals[index] >= 0)
      {
        sample.push_back(_vals[index]);
        _vals[index] = (-_vals[index])-1;
        indices.push_back(index);
        count++;
      }

    } // while(count < sample_size)

    for(int i=0; i<indices.size(); i++)
      _vals[indices[i]] = (-_vals[indices[i]])-1;

  } // select_random_sample(...)

  void select_counter_sample(int sample_size, std::vector<int> &sample)
  {
    int count=0;
   
    while(count < sample_size)
    {
      int index = _distribution(_generator);

      if(_vals[index] >=0)
      {
        _vals[index] = (-_vals[index])-1;
        count++;
      }

    } // while(count < sample_size)

    for(int i=0; i<_vals.size(); i++)
      if(_vals[i]>=0)
        sample.push_back(_vals[i]);
      else
        _vals[i]= (-_vals[i])-1;

  } // select_counter_sample(...)

 private:

  int _n;
  std::vector<int> _vals;
  std::default_random_engine _generator;
  std::uniform_int_distribution<int> _distribution;

}; // class Uniform_sampler

} // PhylogeneticMeasures

#endif // UNIFORM_SAMPLER_H
