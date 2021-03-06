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

#ifndef MEAN_NEAREST_TAXON_DISTANCE_NODE_TYPE_H
#define MEAN_NEAREST_TAXON_DISTANCE_NODE_TYPE_H

#include<iostream>

namespace PhylogeneticMeasures 
{
  template< class KernelType>
  struct Mean_nearest_taxon_distance_auxiliary_data
  {
    typedef KernelType                                          Kernel;
    typedef typename Kernel::Number_type                       Number_type;
    typedef Mean_nearest_taxon_distance_auxiliary_data<Kernel>  Self;

    Number_type first_min, second_min, rest_tree_min;

    Mean_nearest_taxon_distance_auxiliary_data():first_min(Number_type(-1.0)), second_min(Number_type(-1.0)), 
	                                         rest_tree_min(Number_type(-1.0)){}

    Mean_nearest_taxon_distance_auxiliary_data(const Self& d)
    {
      first_min = d.first_min;
      second_min = d.second_min;
      rest_tree_min = d.rest_tree_min;
             
    } // operator=(const Self &d)

    Self& operator=(const Self& d)
    {
      first_min = d.first_min;
      second_min = d.second_min;
      rest_tree_min = d.rest_tree_min;
              
      return *this;
 
    } // operator=(const Self &d)

  }; // struct Mean_nearest_taxon_distance_auxiliary_data


  template< typename KernelType>
  struct Mean_nearest_taxon_distance_node_type:
  public KernelType::template Unimodal_node_augmented<Mean_nearest_taxon_distance_auxiliary_data<KernelType> >
  {};

} // namespace PhylogeneticMeasures 

#endif // MEAN_NEAREST_TAXON_DISTANCE_NODE_TYPE_H
