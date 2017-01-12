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

#ifndef TREE_NODE_BIMODAL_H
#define TREE_NODE_BIMODAL_H

#include<vector>

namespace PhylogeneticMeasures 
{
  template< class KernelType>
  struct Tree_node_bimodal: public KernelType::Unimodal_node
  {
    typedef KernelType                       Kernel;
    typedef typename Kernel::Unimodal_node  Base;
    typedef Tree_node_bimodal<Kernel>        Self;

    bool mark_b; // Auxiliary flag

    int marked_subtree_leaves_b; // Number of species in the subtree
                             // of this node which are also elements
                             // of the query sample.

    std::vector<int> marked_children_b;

    int number_of_marked_children_b()
    { return marked_children_b.size(); }

    Tree_node_bimodal():mark_b(false), marked_subtree_leaves_b(0)
    {}

    Self& operator=(const Self& d)
    {
      Base::operator=(d);

      marked_children_b.clear(); 

      for(int i=0; i<d.marked_children_b.size(); i++)
        marked_children_b.push_back(d.marked_children_b[i]);

      mark_b=d.mark_b; 
      marked_subtree_leaves_b=d.marked_subtree_leaves_b;                

      return *this;
 
   } // operator=(const Self &d)

  }; // struct Tree_node_bimodal

  template< class AUXILIARY_TYPE >
  struct Tree_node_bimodal_augmented: 
  public Tree_node_bimodal<typename AUXILIARY_TYPE::Kernel>, public AUXILIARY_TYPE
  {
    typedef AUXILIARY_TYPE                               Auxiliary_type;
    typedef typename Auxiliary_type::Kernel             Kernel;
    typedef Tree_node_bimodal<Kernel>                    Base_type;
    typedef Tree_node_bimodal_augmented<Auxiliary_type>  Self;

    Tree_node_bimodal_augmented():Base_type(), Auxiliary_type(){}

    Self& operator=(const Self& d)
    {
      Base_type::operator=(d);
      Auxiliary_type::operator=(d);
              
      return *this;
 
   } // operator=(const Self &d)

  }; // Tree_node_bimodal_augmented

} // namespace PhylogeneticMeasures 


#endif //PHYLOGENETIC_TREE_NODE_BIMODAL_H
