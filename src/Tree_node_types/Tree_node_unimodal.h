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

#ifndef TREE_NODE_UNIMODAL_H
#define TREE_NODE_UNIMODAL_H

#include<iostream>
#include<vector>
#include<string>

namespace PhylogeneticMeasures 
{
  // Definitions of a phylogenetic tree node used for measures
  // that involve a single sample of the tree species.

  template <class KernelType>
  struct Tree_node_unimodal
  {
    typedef KernelType                  Kernel;
    typedef Tree_node_unimodal<Kernel>  Self;

    std::string taxon; // Species/taxonomy name.
    double distance; // Distance from father node.
    std::vector<int> children; // Indices to child nodes.
    std::vector<int> marked_children; // Indices to child nodes that have been marked.
    int parent; // Index to parent node.
    bool mark; // Auxiliary flag.
    int all_subtree_leaves; // Number of all leaf nodes.
                        // in the subtree of this node.
    int marked_subtree_leaves; // Number of leaves in the subtree
                           // of this node which are also elements
                           // of the query sample.

    Tree_node_unimodal():distance(-1.0),parent(-1),mark(false),marked_subtree_leaves(0)
    {}

    Tree_node_unimodal(const Self& d)
    {
      this->taxon=d.taxon;
      this->distance=d.distance; 

      this->children.clear();
      this->marked_children.clear(); 

      for(int i=0; i<d.children.size(); i++)
        this->children.push_back(d.children[i]);

      for(int i=0; i<d.marked_children.size(); i++)
        this->marked_children.push_back(d.marked_children[i]);

      this->parent=d.parent; 
      this->mark=d.mark; 
      this->all_subtree_leaves=d.all_subtree_leaves; 
      this->marked_subtree_leaves=d.marked_subtree_leaves;                
 
   } // Tree_node_unimodal(const Self& d)

    int number_of_children() const
    { return children.size(); }

    int number_of_marked_children()
    { return marked_children.size(); }

    Self& operator=(const Self& d)
    {
      taxon=d.taxon;
      distance=d.distance; 

      children.clear();
      marked_children.clear(); 

      for(int i=0; i<d.children.size(); i++)
        children.push_back(d.children[i]);

      for(int i=0; i<d.marked_children.size(); i++)
        marked_children.push_back(d.marked_children[i]);

      parent=d.parent; 
      mark=d.mark; 
      all_subtree_leaves=d.all_subtree_leaves; 
      marked_subtree_leaves=d.marked_subtree_leaves;                

      return *this;
 
   } // operator=(const Self &d)

  }; // struct Tree_node_unimodal


  // Tree node which allows for
  // extra data stored with the node.

  template< class AUXILIARY_TYPE >
  struct Tree_node_unimodal_augmented: 
  public Tree_node_unimodal<typename AUXILIARY_TYPE::Kernel>, public AUXILIARY_TYPE
  {
    typedef AUXILIARY_TYPE                                Auxiliary_type;
    typedef typename Auxiliary_type::Kernel              Kernel;
    typedef Tree_node_unimodal<Kernel>                    Base_type;
    typedef Tree_node_unimodal_augmented<Auxiliary_type>  Self;

    Tree_node_unimodal_augmented():Base_type(), Auxiliary_type(){}

    Tree_node_unimodal_augmented(const Self &d):Base_type(d),Auxiliary_type(d){}


    Self& operator=(const Self& d)
    {
      Base_type::operator=(d);
      Auxiliary_type::operator=(d);
              
      return *this;
 
   } // operator=(const Self &d)

  }; // Tree_node_unimodal_augmented

} // namespace PhylogeneticMeasures 

#endif //PHYLOGENETIC_TREE_NODE_UNIMODAL_H
