//////////////////////////////////////////////////////////////////////////////////
//    Copyright (C) 2015,  Constantinos Tsirogiannis.  Email: analekta@gmail.com
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

#ifndef EXCEPTION_FUNCTOR_SILENT_H
#define EXCEPTION_FUNCTOR_SILENT_H

#include<string>
#include<cstdlib>
#include<exception>
#include<vector>
#include"Exception_type.h"


class Warning_list_type
{
 public: 

  int number_of_warnings()
  { return _warnings.size();}

  std::string operator[](int i)
  { return _warnings[i];}

  void push_warning(std::string &w)
  { _warnings.push_back(w);} 

  void clear()
  { _warnings.clear();}

 private:

  std::vector<std::string> _warnings;
}; 

Warning_list_type warning_list;

namespace ExceptionRelatedTypes
{
  class Exception_functor
  {
   public:

    Exception_functor(){}
    
    void operator()(Exception_type excp)
    { throw excp;}

    void issue_warning(std::string str)
    { warning_list.push_warning(str);}

  }; // class Exception_functor

} // ExceptionRelatedTypes

#endif // EXCEPTION_FUNCTOR_SILENT_H
