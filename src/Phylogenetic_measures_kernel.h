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

#ifndef PHYLOGENETIC_MEASURES_KERNEL_H
#define PHYLOGENETIC_MEASURES_KERNEL_H

#include "Tree_node.h"
#include "Measures/Measure_base.h"
#include "Measures/Polynomial_related_types.h"
#include "Measures/Poisson_binomial_types.h"
#include "Incremental_Monte_Carlo_types/Incremental_Monte_Carlo_handler.h"
#include "Random_samplers.h"
#include "Phylogenetic_measures.h"
#include "Numeric_traits_double.h"
#include "Phylogenetic_tree_base.h"
#include "Phylogenetic_tree_bimodal.h"
#include "Exception_related_types.h"

template< class NTS = typename PhylogeneticMeasures::Numeric_traits_double>
struct Phylogenetic_measures_kernel
{
  typedef NTS                                              Numeric_traits;
  typedef Phylogenetic_measures_kernel<Numeric_traits>     Self;
  typedef typename Numeric_traits::Number_type            Number_type;
  typedef typename Numeric_traits::Protected_number_type  Protected_number_type;
  
  enum Measure_type {PHYLOGENETIC_DIVERSITY, MEAN_PAIRWISE_DISTANCE, 
                     MEAN_NEAREST_TAXON_DISTANCE, CORE_ANCESTOR_COST, 
                     COMMON_BRANCH_LENGTH, COMMUNITY_DISTANCE, 
                     COMMUNITY_DISTANCE_NEAREST_TAXON, 
                     PHYLOGENETIC_SORENSENS_SIMILARITY, 
                     UNIQUE_FRACTION};

  enum Distribution_type {UNIFORM_FIXED_SIZE, POISSON_BINOMIAL, POISSON_BINOMIAL_FIXED_SIZE, SEQUENTIAL_FIXED_SIZE};
  enum Edge_relation_type {OFFSPRING, ANCESTOR, INDEPENDENT};
  
  /////////////////////////////////////////////////////////////////////////////
  // Node types that can be used for implementing a phylogenetic tree class. //
  /////////////////////////////////////////////////////////////////////////////
  
  typedef typename PhylogeneticMeasures::Tree_node_unimodal<Self>      Unimodal_node;
  typedef typename PhylogeneticMeasures::Tree_node_bimodal<Self>       Bimodal_node;
  
  template <class AUXILIARY_TYPE> 
  struct Unimodal_node_augmented:
  public PhylogeneticMeasures::Tree_node_unimodal_augmented<AUXILIARY_TYPE>
  {};

  template <class AUXILIARY_TYPE> 
  struct Bimodal_node_augmented:
  public PhylogeneticMeasures::Tree_node_bimodal_augmented<AUXILIARY_TYPE>
  {};

  typedef typename PhylogeneticMeasures::Mean_nearest_taxon_distance_node_type<Self>  
                                                                 Mean_nearest_taxon_distance_node_type;
  
  typedef typename PhylogeneticMeasures::Community_distance_nearest_taxon_node_type<Self>  
                                                                 Community_distance_nearest_taxon_node_type;  
   
  ////////////////////////////////
  // Phylogenetic tree classes  //
  ////////////////////////////////
  
  typedef typename PhylogeneticMeasures::Phylogenetic_tree_base<Self,Unimodal_node>    Unimodal_tree;
  typedef typename PhylogeneticMeasures::Phylogenetic_tree_bimodal<Self,Bimodal_node>  Bimodal_tree;
  
  typedef typename PhylogeneticMeasures::Phylogenetic_tree_base<Self,Mean_nearest_taxon_distance_node_type>
                                                                                Mean_nearest_taxon_distance_tree;
										 
  typedef typename PhylogeneticMeasures::Phylogenetic_tree_bimodal<Self,Community_distance_nearest_taxon_node_type>
                                                                                Community_distance_nearest_taxon_tree;  

																				
  ////////////////////////////////////////////////////////////////////////////////
  // Base types for implementing classes of phylogenetic biodiversity measures  //
  ////////////////////////////////////////////////////////////////////////////////
																				
  typedef typename PhylogeneticMeasures::Measure_base_unimodal<Self>  Measure_base_unimodal;
  typedef typename PhylogeneticMeasures::Measure_base_bimodal<Self>   Measure_base_bimodal;

  template <class TreeType> 
  struct Mean_pairwise_distance_base:
  public PhylogeneticMeasures::Mean_pairwise_distance_base<Self,TreeType>
  {};

  ///////////////////////////////////////////////////////////////////
  // Classes that are used to compute the moments of several       // 
  // measures under the constrained Poisson binomial distribution  //
  ///////////////////////////////////////////////////////////////////

  typedef typename PhylogeneticMeasures::Poisson_binomial_moments_Phylogenetic_diversity<Self>
                                                  Poisson_binomial_moments_Phylogenetic_diversity;

  typedef typename PhylogeneticMeasures::Poisson_binomial_moments_Mean_pairwise_distance<Self>
                                                  Poisson_binomial_moments_Mean_pairwise_distance;

  typedef typename PhylogeneticMeasures::Poisson_binomial_moments_Mean_nearest_taxon_distance<Self>
                                                  Poisson_binomial_moments_Mean_nearest_taxon_distance;

  ////////////////////////////////////////////////////////////////////////////////////////
  // Classes that compute the values and moments of phylogenetic biodiversity measures  //
  ////////////////////////////////////////////////////////////////////////////////////////
  
  typedef typename PhylogeneticMeasures::Phylogenetic_diversity<Self>             Phylogenetic_diversity;
  typedef typename PhylogeneticMeasures::Mean_pairwise_distance<Self>             Mean_pairwise_distance;
  typedef typename PhylogeneticMeasures::Mean_nearest_taxon_distance<Self>        Mean_nearest_taxon_distance;
  typedef typename PhylogeneticMeasures::Core_ancestor_cost<Self>                 Core_ancestor_cost;
  typedef typename PhylogeneticMeasures::Common_branch_length<Self>               Common_branch_length;
  typedef typename PhylogeneticMeasures::Community_distance<Self>                 Community_distance;
  typedef typename PhylogeneticMeasures::Community_distance_nearest_taxon<Self>   Community_distance_nearest_taxon;  
  typedef typename PhylogeneticMeasures::Phylogenetic_Sorensens_similarity<Self>  Phylogenetic_Sorensens_similarity;
  typedef typename PhylogeneticMeasures::Unique_fraction<Self>                    Unique_fraction;


  ///////////////////////////////////////////////////
  // Classes related to polynomial multiplication  //
  ///////////////////////////////////////////////////

  typedef typename PhylogeneticMeasures::Polynomial_rep<Self>             Polynomial;
  typedef typename PhylogeneticMeasures::Polynomial_multiplication<Self>  Polynomial_multiplication;

  ///////////////////////////////////////////
  // Classes that perform random sampling. //
  ///////////////////////////////////////////

  typedef typename PhylogeneticMeasures::Uniform_sampler<Self>     Uniform_sampler;
  typedef typename PhylogeneticMeasures::Sequential_sampler<Self>  Sequential_sampler;

  ///////////////////////////////////////////////////////////////
  // Class that performs incremental Monte Carlo calculations. //
  ///////////////////////////////////////////////////////////////

  typedef typename PhylogeneticMeasures::Incremental_Monte_Carlo_handler<Self> 
                                                              Incremental_Monte_Carlo_handler;

  //////////////////////////////////////////
  // Classes used for handling exceptions //
  //////////////////////////////////////////

  typedef typename ExceptionRelatedTypes::Exception_functor  Exception_functor;
  typedef typename ExceptionRelatedTypes::Exception_type     Exception_type;  

}; // class Phylogenetic_measures_kernel

#endif //PHYLOGENETIC_MEASURES_KERNEL_H
