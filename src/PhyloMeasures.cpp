////////////////////////////////////////////////////////////////////////////////////////////////
//    Copyright (C) 2016,  Constantinos Tsirogiannis and Brody Sandel.
//
//    Email: tsirogiannis.c@gmail.com and bsandel@scu.edu
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
////////////////////////////////////////////////////////////////////////////////////////////////

#include<R.h>
#include<R_ext/Print.h>
#include<vector>
#include<iostream>
#include"Phylogenetic_measures_kernel.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// General Types ////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef PhylogeneticMeasures::Numeric_traits_double                     Numeric_traits;  
typedef Phylogenetic_measures_kernel<Numeric_traits>                    Kernel;
typedef Kernel::Number_type                                             Number_type;
typedef Numeric_traits::Square_root                                     Square_root;

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Tree Types /////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef Kernel::Unimodal_tree                                           Unimodal_tree;
typedef Kernel::Bimodal_tree                                            Bimodal_tree;
typedef Kernel::Mean_nearest_taxon_distance_tree                        Mean_nearest_taxon_distance_tree;
typedef Kernel::Community_distance_nearest_taxon_tree                   Community_distance_nearest_taxon_tree;

///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// Measures Types ////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef Kernel::Mean_pairwise_distance                                  Mean_pairwise_distance;
typedef Kernel::Phylogenetic_diversity                                  Phylogenetic_diversity;
typedef Kernel::Core_ancestor_cost                                      Core_ancestor_cost;
typedef Kernel::Mean_nearest_taxon_distance                             Mean_nearest_taxon_distance;
typedef Kernel::Common_branch_length                                    Common_branch_length;
typedef Kernel::Community_distance                                      Community_distance;
typedef Kernel::Community_distance_nearest_taxon                        Community_distance_nearest_taxon;
typedef Kernel::Phylogenetic_Sorensens_similarity                       Phylogenetic_Sorensens_similarity;
typedef Kernel::Unique_fraction                                         Unique_fraction;

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////// Moments and p-Value Computations Involving Abundance Weighting //////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef Kernel::Poisson_binomial_moments_Phylogenetic_diversity         Poisson_PD;
typedef Kernel::Poisson_binomial_moments_Mean_pairwise_distance         Poisson_MPD;
typedef Kernel::Poisson_binomial_moments_Mean_nearest_taxon_distance    Poisson_MNTD;

typedef Kernel::Sequential_sampler                                      Sequential_sampler;
typedef Kernel::Incremental_Monte_Carlo_handler                         Incremental_Monte_Carlo_handler;

///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Exception Handling Types //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef Kernel::Exception_type                                          Exception_type;

extern Warning_list_type warning_list;

extern "C" {


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///////////////// Helper functions /////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

void flush_warnings()
{
  if(warning_list.number_of_warnings() > 0)
  {
    REprintf("\n");

    for(int i=0; i<warning_list.number_of_warnings(); i++)
      REprintf("%s \n", warning_list[i].c_str()); 

    REprintf("\n");

    warning_list.clear();
  }
}

void copy_message(std::string &msg, char **message)
{
  for(int i=0; i<msg.size(); i++)
    (*message)[i] = msg[i];
}

void transform_matrix_query_arguments_unimodal( int *n_w, int *n_l, int *edges, double *weights, 
                                                char **names, char **matrix_names, int *n_r, int *n_c,
                                                int *matrix_data, std::vector<int> &parents, 
                                                std::vector<int> &children, std::vector<Number_type> &edge_weights,
                                                std::vector<std::string> &species_names,
                                                std::vector<std::string> &matrix_species_names, 
                                                std::vector<std::vector<bool> > &matrix )
{
  int number_of_nodes = (*n_w)+1,
      number_of_edges = (*n_w), 
      number_of_leaves = *n_l;

  for(int i=0; i<(*n_w); i++)
    parents.push_back(edges[i]);

  for(int i=(*n_w); i<2*(*n_w); i++)
    children.push_back(edges[i]);

  for(int i=0; i<number_of_edges; i++)
    edge_weights.push_back(weights[i]);

  for(int i=0; i<number_of_leaves; i++)
    species_names.push_back(names[i]);

  for(int i=0; i<*n_c; i++)
    matrix_species_names.push_back(matrix_names[i]);

  matrix.assign(*n_r, std::vector<bool>()); 

  for(int i=0; i<(*n_r)*(*n_c); i++)
  {
    int row = i/(*n_c);
 
    matrix[row].push_back(bool(matrix_data[i]));
  }  

} // transform_matrix_query_arguments_unimodal(...)


void transform_matrix_query_arguments_bimodal( int *n_w, int *n_l, int *edges, double *weights, 
                                               char **names, char **matrix_names_a, int *n_r_a, int *n_c_a, 
                                               int *matrix_data_a, char **matrix_names_b, int *n_r_b, int *n_c_b,  
                                               int *matrix_data_b, int *number_of_pairs, int *pair_data, 
                                               std::vector<int> &parents, std::vector<int> &children,
                                               std::vector<Number_type> &edge_weights,
                                               std::vector<std::string> &species_names,
                                               std::vector<std::string> &matrix_species_names_a, 
                                               std::vector<std::vector<bool> > &matrix_a,
                                               std::vector<std::string> &matrix_species_names_b, 
                                               std::vector<std::vector<bool> > &matrix_b,
                                               std::vector<std::pair<int, int> > &query_pairs)
{
  int number_of_nodes = (*n_w)+1,
      number_of_edges = (*n_w), 
      number_of_leaves = *n_l;

  for(int i=0; i<(*n_w); i++)
    parents.push_back(edges[i]);

  for(int i=(*n_w); i<2*(*n_w); i++)
    children.push_back(edges[i]);

  for(int i=0; i<number_of_edges; i++)
    edge_weights.push_back(weights[i]);

  for(int i=0; i<number_of_leaves; i++)
    species_names.push_back(names[i]);

  for(int i=0; i<*n_c_a; i++)
    matrix_species_names_a.push_back(matrix_names_a[i]);

  matrix_a.assign(*n_r_a, std::vector<bool>());

  for(int i=0; i<(*n_r_a)*(*n_c_a); i++)
  {
    int row = i/(*n_c_a);
 
    matrix_a[row].push_back(bool(matrix_data_a[i]));
  }  

  if(*n_r_b > 0)
  {
    for(int i=0; i<*n_c_b; i++)
      matrix_species_names_b.push_back(matrix_names_b[i]);

    matrix_b.assign(*n_r_b, std::vector<bool>()); 

    for(int i=0; i<(*n_r_b)*(*n_c_b); i++)
    {
      int row = i/(*n_c_b);
 
      matrix_b[row].push_back(bool(matrix_data_b[i]));
    }
  } 

  if(*number_of_pairs > 0)
  {
    query_pairs.assign(*number_of_pairs, std::make_pair(-1,-1));

    for(int i=0; i<2*(*number_of_pairs); i+=2)
    {
      query_pairs[i/2].first=pair_data[i]-1;
      query_pairs[i/2].second=pair_data[i+1]-1;
    }
  }

} // transform_matrix_query_arguments_bimodal(...)

void transform_moments_function_arguments_unimodal( int *n_w, int *n_l, int *edges, double *weights, 
                                                    char **names, int *n_q, int *q_sizes, 
                                                    std::vector<int> &parents, std::vector<int> &children,
                                                    std::vector<Number_type> &edge_weights,
                                                    std::vector<std::string> &species_names,
                                                    std::vector<int> &query_sizes)
{
  int number_of_nodes = (*n_w)+1,
      number_of_edges = (*n_w), 
      number_of_leaves = *n_l;

  for(int i=0; i<(*n_w); i++)
    parents.push_back(edges[i]);

  for(int i=(*n_w); i<2*(*n_w); i++)
    children.push_back(edges[i]);

  for(int i=0; i<number_of_edges; i++)
    edge_weights.push_back(weights[i]);

  for(int i=0; i<number_of_leaves; i++)
    species_names.push_back(names[i]);

  for(int i=0; i<*n_q; i++)
    query_sizes.push_back(q_sizes[i]);

} // transform_moments_function_arguments_unimodal(...)

void transform_moments_function_arguments_bimodal( int *n_w, int *n_l, int *edges, double *weights, 
                                                   char **names, int *n_q, int *q_sizes, 
                                                   std::vector<int> &parents, std::vector<int> &children,
                                                   std::vector<Number_type> &edge_weights,
                                                   std::vector<std::string> &species_names,
                                                   std::vector<std::pair<int,int> > &query_sizes)
{
  int number_of_nodes = (*n_w)+1,
      number_of_edges = (*n_w), 
      number_of_leaves = *n_l;

  for(int i=0; i<(*n_w); i++)
    parents.push_back(edges[i]);

  for(int i=(*n_w); i<2*(*n_w); i++)
    children.push_back(edges[i]);

  for(int i=0; i<number_of_edges; i++)
    edge_weights.push_back(weights[i]);

  for(int i=0; i<number_of_leaves; i++)
    species_names.push_back(names[i]);

  int actual_size = *n_q/2;

  for(int i=0; i<*n_q/2; i++)
  {
    std::pair<int,int> temp;
    temp.first = q_sizes[i];
    temp.second = q_sizes[actual_size+i];
    query_sizes.push_back(temp);
  }

} // transform_moments_function_arguments_bimodal(...)

void transform_abundance_weights( int *n_l, char **abundance_weights_nms, 
                                  double *abundance_weights_vals, 
                                  std::vector<std::string> &abundance_weights_names,
                                  std::vector<Number_type> &abundance_weights_values)
{
  int number_of_leaves = *n_l;

  for(int i=0; i<number_of_leaves; i++)
    abundance_weights_values.push_back(abundance_weights_vals[i]);

  for(int i=0; i<number_of_leaves; i++)
    abundance_weights_names.push_back(abundance_weights_nms[i]);

} // transform_abundance_weights(...)


////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////// Matrix query functions //////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

void pd_query(int *n_w, int *n_l, int *edges, double *weights, char **names,
              char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
              bool *standardised, double *output,
              char **error_message, int *message_size)
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Phylogenetic_diversity pd(tree);

    std::vector<Number_type> result;

    if(*standardised == true)
      pd.matrix_query_standardised(matrix_species_names, matrix, std::back_inserter(result));  
    else
      pd.matrix_query_basic(matrix_species_names, matrix, std::back_inserter(result));    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void pd_query(...)


void mpd_query(int *n_w, int *n_l, int *edges, double *weights, char **names,
               char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
               bool *standardised, double *output,
               char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Mean_pairwise_distance mpd(tree);

    std::vector<Number_type> result;

    if(*standardised == true)
      mpd.matrix_query_standardised(matrix_species_names, matrix, std::back_inserter(result));  
    else
      mpd.matrix_query_basic(matrix_species_names, matrix, std::back_inserter(result));    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void mpd_query(...)


void mntd_query(int *n_w, int *n_l, int *edges, double *weights, char **names,
                char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                bool *standardised, double *output,
                char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  try
  {
    Mean_nearest_taxon_distance_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Mean_nearest_taxon_distance mntd(tree);

    std::vector<Number_type> result;

    if(*standardised == true)
      mntd.matrix_query_standardised(matrix_species_names, matrix, std::back_inserter(result));  
    else
      mntd.matrix_query_basic(matrix_species_names, matrix, std::back_inserter(result));    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void mntd_query(...)


/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

void cac_query(int *n_w, int *n_l, int *edges, double *weights, char **names, double *chi,
               char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
               bool *standardised, double *output,
               char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Core_ancestor_cost cac(tree, *chi);

    std::vector<Number_type> result;

    if(*standardised == true)
      cac.matrix_query_standardised(matrix_species_names, matrix, std::back_inserter(result));  
    else
      cac.matrix_query_basic(matrix_species_names, matrix, std::back_inserter(result));    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void cac_query(...)

////////////////////////////////////////////////////
////////////////////////////////////////////////////
//// Abundance weighted matrix query functions /////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

void pd_query_abundance_weighted(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                 char **abundance_weights_nms, double *abundance_weights_vals, 
                                 char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                                 bool *standardised, double *output,
                                 char **error_message, int *message_size)
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights,abundance_weights_values;
  std::vector<std::string> species_names, matrix_species_names, abundance_weights_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Phylogenetic_diversity pd(tree);

    std::vector<Number_type> result;
 
    pd.set_probability_distribution(Kernel::POISSON_BINOMIAL_FIXED_SIZE);

    if(*standardised == true)
      pd.matrix_query_weighted_standardised(matrix_species_names, matrix, std::back_inserter(result));  
    else
      pd.matrix_query_weighted_basic(matrix_species_names, matrix, std::back_inserter(result));    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void pd_query_abundance_weighted(...)


void mpd_query_abundance_weighted(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                  char **abundance_weights_nms, double *abundance_weights_vals, 
                                  char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                                  bool *standardised, double *output,
                                  char **error_message, int *message_size)
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights,abundance_weights_values;
  std::vector<std::string> species_names, matrix_species_names, abundance_weights_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Mean_pairwise_distance mpd(tree);

    std::vector<Number_type> result;
 
    mpd.set_probability_distribution(Kernel::POISSON_BINOMIAL_FIXED_SIZE);

    if(*standardised == true)
      mpd.matrix_query_weighted_standardised(matrix_species_names, matrix, std::back_inserter(result));  
    else
      mpd.matrix_query_weighted_basic(matrix_species_names, matrix, std::back_inserter(result));    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void mpd_query_abundance_weighted(...)

void mntd_query_abundance_weighted(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                   char **abundance_weights_nms, double *abundance_weights_vals, 
                                   char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                                   bool *standardised, double *output,
                                   char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights, abundance_weights_values;
  std::vector<std::string> species_names, matrix_species_names,abundance_weights_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Mean_nearest_taxon_distance_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Mean_nearest_taxon_distance mntd(tree);

    mntd.set_probability_distribution(Kernel::POISSON_BINOMIAL_FIXED_SIZE);

    std::vector<Number_type> result;

    if(*standardised == true)
      mntd.matrix_query_weighted_standardised(matrix_species_names, matrix, std::back_inserter(result));  
    else
      mntd.matrix_query_weighted_basic(matrix_species_names, matrix, std::back_inserter(result));    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void mntd_query_abundance_weighted(...)


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//// Abundance-weighted matrix query functions that use sequential Monte-Carlo computations /////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void pd_query_weighted_sequential(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                  char **abundance_weights_nms, double *abundance_weights_vals, 
                                  char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                                  bool *standardised, int *repetitions, int *seed, double *output,
                                  char **error_message, int *message_size)
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights,abundance_weights_values;
  std::vector<std::string> species_names, matrix_species_names, abundance_weights_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Phylogenetic_diversity pd(tree);

    std::vector<Number_type> result;
    int reps = *repetitions; 

    pd.set_probability_distribution(Kernel::SEQUENTIAL_FIXED_SIZE);
    pd.set_seed(*seed);

    if(*standardised == true)
      pd.matrix_query_weighted_standardised(matrix_species_names, matrix, 
                                            std::back_inserter(result), reps);  
    else
      pd.matrix_query_weighted_basic(matrix_species_names, matrix, std::back_inserter(result));    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void pd_query_weighted_sequential(...)


void mpd_query_weighted_sequential(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                   char **abundance_weights_nms, double *abundance_weights_vals, 
                                   char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                                   bool *standardised, int *repetitions, int *seed, double *output,
                                   char **error_message, int *message_size)
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights,abundance_weights_values;
  std::vector<std::string> species_names, matrix_species_names, abundance_weights_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Mean_pairwise_distance mpd(tree);

    std::vector<Number_type> result;
    int reps = *repetitions; 

    mpd.set_probability_distribution(Kernel::SEQUENTIAL_FIXED_SIZE);
    mpd.set_seed(*seed);

    if(*standardised == true)
      mpd.matrix_query_weighted_standardised(matrix_species_names, matrix, 
                                             std::back_inserter(result), reps);  
    else
      mpd.matrix_query_weighted_basic(matrix_species_names, matrix, std::back_inserter(result));    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void mpd_query_weighted_sequential(...)

void mntd_query_weighted_sequential(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                    char **abundance_weights_nms, double *abundance_weights_vals, 
                                    char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                                    bool *standardised, int *repetitions, int *seed, double *output,
                                    char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights, abundance_weights_values;
  std::vector<std::string> species_names, matrix_species_names,abundance_weights_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Mean_nearest_taxon_distance_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Mean_nearest_taxon_distance mntd(tree);

    mntd.set_probability_distribution(Kernel::SEQUENTIAL_FIXED_SIZE);
    mntd.set_seed(*seed);

    std::vector<Number_type> result;
    int reps = *repetitions;

    if(*standardised == true)
      mntd.matrix_query_weighted_standardised(matrix_species_names, matrix, 
                                              std::back_inserter(result),reps);  
    else
      mntd.matrix_query_weighted_basic(matrix_species_names, matrix, std::back_inserter(result));    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void mntd_query_weighted_sequential(...)


void cac_query_weighted_sequential(int *n_w, int *n_l, int *edges, double *weights,  
                                    char **names, double *chi, char **abundance_weights_nms, 
                                    double *abundance_weights_vals, char **matrix_names, 
                                    int *n_r, int *n_c, int *matrix_data, 
                                    bool *standardised, int *repetitions, int *seed, double *output,
                                    char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights, abundance_weights_values;
  std::vector<std::string> species_names, matrix_species_names,abundance_weights_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Core_ancestor_cost cac(tree,*chi);

    cac.set_probability_distribution(Kernel::SEQUENTIAL_FIXED_SIZE);
    cac.set_seed(*seed);

    std::vector<Number_type> result;
    int reps = *repetitions;

    if(*standardised == true)
      cac.matrix_query_weighted_standardised(matrix_species_names, matrix, 
                                              std::back_inserter(result),reps);  
    else
      cac.matrix_query_weighted_basic(matrix_species_names, matrix, std::back_inserter(result));    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void cac_query_weighted_sequential(...)


void cbl_query(int *n_w, int *n_l, int *edges, double *weights, char **names,
               char **matrix_names_a, int *n_r_a, int *n_c_a, int *matrix_data_a, 
               char **matrix_names_b, int *n_r_b, int *n_c_b, int *matrix_data_b,
               int *n_p, int *pair_data, bool *standardised, double *output,
               char **error_message, int *message_size)
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names_a, matrix_species_names_b;
  std::vector<std::vector<bool> > matrix_a, matrix_b;
  std::vector<std::pair<int,int> > query_pairs;  

  transform_matrix_query_arguments_bimodal(n_w, n_l, edges, weights, names, 
                                           matrix_names_a, n_r_a, n_c_a, matrix_data_a,
                                           matrix_names_b, n_r_b, n_c_b, matrix_data_b,
                                           n_p, pair_data, parents, children, 
                                           edge_weights, species_names, 
                                           matrix_species_names_a, matrix_a,
                                           matrix_species_names_b, matrix_b, 
                                           query_pairs);
  try
  {
    Bimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Common_branch_length cbl(tree);

    std::vector<Number_type> result;

    if(matrix_b.size() > 0)
    {
      if(query_pairs.size() > 0)
      {
        if(*standardised == true)
          cbl.matrix_query_specific_pairs_standardised(matrix_species_names_a, matrix_a,
                                                       matrix_species_names_b, matrix_b, 
                                                       query_pairs, std::back_inserter(result));  
        else
          cbl.matrix_query_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                                matrix_species_names_b, matrix_b, 
                                                query_pairs, std::back_inserter(result)); 
      } // if(query_pairs.size() > 0)
      else
      {
        if(*standardised == true)
          cbl.matrix_query_standardised(matrix_species_names_a, matrix_a,
                                        matrix_species_names_b, matrix_b, 
                                        std::back_inserter(result));  
        else
          cbl.matrix_query_basic(matrix_species_names_a, matrix_a,
                                 matrix_species_names_b, matrix_b, 
                                 std::back_inserter(result));

      } // else of if(query_pairs.size() > 0)   

    } // if(matrix_b.size() > 0)
    else
    {
      if(query_pairs.size() > 0)
      {
        if(*standardised == true)
          cbl.matrix_query_specific_pairs_standardised(matrix_species_names_a, matrix_a,
                                                       query_pairs, std::back_inserter(result));  
        else
          cbl.matrix_query_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                                query_pairs, std::back_inserter(result)); 
      } // if(query_pairs.size() > 0)
      else
      {
        if(*standardised == true)
          cbl.matrix_query_standardised(matrix_species_names_a, matrix_a,
                                        std::back_inserter(result));  
        else
          cbl.matrix_query_basic(matrix_species_names_a, matrix_a,
                                 std::back_inserter(result));

      } // else of if(query_pairs.size() > 0)   

    } // else of if(matrix_b.size() > 0) 

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void cbl_query(...)


void cd_query(int *n_w, int *n_l, int *edges, double *weights, char **names,
              char **matrix_names_a, int *n_r_a, int *n_c_a, int *matrix_data_a, 
              char **matrix_names_b, int *n_r_b, int *n_c_b, int *matrix_data_b,
              int *n_p, int *pair_data, bool *standardised, double *output,
              char **error_message, int *message_size)
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names_a, matrix_species_names_b;
  std::vector<std::vector<bool> > matrix_a, matrix_b;
  std::vector<std::pair<int,int> > query_pairs;  

  transform_matrix_query_arguments_bimodal(n_w, n_l, edges, weights, names, 
                                           matrix_names_a, n_r_a, n_c_a, matrix_data_a,
                                           matrix_names_b, n_r_b, n_c_b, matrix_data_b,
                                           n_p, pair_data, parents, children, 
                                           edge_weights, species_names, 
                                           matrix_species_names_a, matrix_a,
                                           matrix_species_names_b, matrix_b, 
                                           query_pairs);

  try
  {
    Bimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Community_distance cd(tree);

    std::vector<Number_type> result;

    if(matrix_b.size() > 0)
    {
      if(query_pairs.size() > 0)
      {
        if(*standardised == true)
          cd.matrix_query_specific_pairs_standardised(matrix_species_names_a, matrix_a,
                                                      matrix_species_names_b, matrix_b, 
                                                      query_pairs, std::back_inserter(result));  
        else
          cd.matrix_query_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                               matrix_species_names_b, matrix_b, 
                                               query_pairs, std::back_inserter(result)); 
      } // if(query_pairs.size() > 0)
      else
      {
        if(*standardised == true)
          cd.matrix_query_standardised(matrix_species_names_a, matrix_a,
                                       matrix_species_names_b, matrix_b, 
                                       std::back_inserter(result));  
        else
          cd.matrix_query_basic(matrix_species_names_a, matrix_a,
                                matrix_species_names_b, matrix_b, 
                                std::back_inserter(result));

      } // else of if(query_pairs.size() > 0)   

    } // if(matrix_b.size() > 0)
    else
    {
      if(query_pairs.size() > 0)
      {
        if(*standardised == true)
          cd.matrix_query_specific_pairs_standardised(matrix_species_names_a, matrix_a,
                                                      query_pairs, std::back_inserter(result));  
        else
          cd.matrix_query_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                               query_pairs, std::back_inserter(result)); 
      } // if(query_pairs.size() > 0)
      else
      {
        if(*standardised == true)
          cd.matrix_query_standardised(matrix_species_names_a, matrix_a,
                                       std::back_inserter(result));  
        else
          cd.matrix_query_basic(matrix_species_names_a, matrix_a,
                                std::back_inserter(result));

      } // else of if(query_pairs.size() > 0)   

    } // else of if(matrix_b.size() > 0) 


    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void cd_query(...)


void cdnt_query(int *n_w, int *n_l, int *edges, double *weights, char **names,
                char **matrix_names_a, int *n_r_a, int *n_c_a, int *matrix_data_a, 
                char **matrix_names_b, int *n_r_b, int *n_c_b, int *matrix_data_b,
                int *n_p, int *pair_data, double *output,
                char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names_a, matrix_species_names_b;
  std::vector<std::vector<bool> > matrix_a, matrix_b;
  std::vector<std::pair<int,int> > query_pairs;  

  transform_matrix_query_arguments_bimodal(n_w, n_l, edges, weights, names, 
                                           matrix_names_a, n_r_a, n_c_a, matrix_data_a,
                                           matrix_names_b, n_r_b, n_c_b, matrix_data_b,
                                           n_p, pair_data, parents, children, 
                                           edge_weights, species_names, 
                                           matrix_species_names_a, matrix_a,
                                           matrix_species_names_b, matrix_b, 
                                           query_pairs);

  try
  {
    Community_distance_nearest_taxon_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Community_distance_nearest_taxon cdnt(tree);
 
    std::vector<Number_type> result;

    if(matrix_b.size() > 0)
    {
      if(query_pairs.size() > 0)
      {
        cdnt.matrix_query_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                               matrix_species_names_b, matrix_b, 
                                               query_pairs, std::back_inserter(result)); 
      }
      else
      {
        cdnt.matrix_query_basic(matrix_species_names_a, matrix_a,
                                matrix_species_names_b, matrix_b, 
                                std::back_inserter(result));
      }   

    } // if(matrix_b.size() > 0)
    else
    {
      if(query_pairs.size() > 0)
      {
        cdnt.matrix_query_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                               query_pairs, std::back_inserter(result)); 
      } 
      else
      {
        cdnt.matrix_query_basic(matrix_species_names_a, matrix_a,
                                std::back_inserter(result));
      }   

    } // else of if(matrix_b.size() > 0) 


    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void cdnt_query(...)


void cdnt_averaged_query(int *n_w, int *n_l, int *edges, double *weights, char **names,
                         char **matrix_names_a, int *n_r_a, int *n_c_a, int *matrix_data_a, 
                         char **matrix_names_b, int *n_r_b, int *n_c_b, int *matrix_data_b,
                         int *n_p, int *pair_data, double *output,
                         char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names_a, matrix_species_names_b;
  std::vector<std::vector<bool> > matrix_a, matrix_b;
  std::vector<std::pair<int,int> > query_pairs;  

  transform_matrix_query_arguments_bimodal(n_w, n_l, edges, weights, names, 
                                           matrix_names_a, n_r_a, n_c_a, matrix_data_a,
                                           matrix_names_b, n_r_b, n_c_b, matrix_data_b,
                                           n_p, pair_data, parents, children, 
                                           edge_weights, species_names, 
                                           matrix_species_names_a, matrix_a,
                                           matrix_species_names_b, matrix_b, 
                                           query_pairs);

  try
  {
    Community_distance_nearest_taxon_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Community_distance_nearest_taxon cdnt(tree);

    std::vector<Number_type> result;

    if(matrix_b.size() > 0)
    {
      if(query_pairs.size() > 0)
      {
        cdnt.matrix_query_averaged_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                                        matrix_species_names_b, matrix_b, 
                                                        query_pairs, std::back_inserter(result)); 
      }
      else
      {
        cdnt.matrix_query_averaged_basic(matrix_species_names_a, matrix_a,
                                         matrix_species_names_b, matrix_b, 
                                         std::back_inserter(result));
      }   

    } // if(matrix_b.size() > 0)
    else
    {
      if(query_pairs.size() > 0)
      {
        cdnt.matrix_query_averaged_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                                        query_pairs, std::back_inserter(result)); 
      } 
      else
      {
        cdnt.matrix_query_averaged_basic(matrix_species_names_a, matrix_a,
                                         std::back_inserter(result));
      }   

    } // else of if(matrix_b.size() > 0) 


    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void cdnt_query_averaged(...)


void cdnt_directed_query(int *n_w, int *n_l, int *edges, double *weights, char **names,
                         char **matrix_names_a, int *n_r_a, int *n_c_a, int *matrix_data_a, 
                         char **matrix_names_b, int *n_r_b, int *n_c_b, int *matrix_data_b,
                         int *n_p, int *pair_data, double *output,
                         char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names_a, matrix_species_names_b;
  std::vector<std::vector<bool> > matrix_a, matrix_b;
  std::vector<std::pair<int,int> > query_pairs;  

  transform_matrix_query_arguments_bimodal(n_w, n_l, edges, weights, names, 
                                           matrix_names_a, n_r_a, n_c_a, matrix_data_a,
                                           matrix_names_b, n_r_b, n_c_b, matrix_data_b,
                                           n_p, pair_data, parents, children, 
                                           edge_weights, species_names, 
                                           matrix_species_names_a, matrix_a,
                                           matrix_species_names_b, matrix_b, 
                                           query_pairs);

  try
  {
    Community_distance_nearest_taxon_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Community_distance_nearest_taxon cdnt(tree);

    std::vector<Number_type> result_a_to_b, result_b_to_a;

    if(matrix_b.size() > 0)
    {
      if(query_pairs.size() > 0)
      {
        cdnt.matrix_query_directed_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                                        matrix_species_names_b, matrix_b, 
                                                        query_pairs, std::back_inserter(result_a_to_b),
                                                        std::back_inserter(result_b_to_a)); 
      }
      else
      {
        cdnt.matrix_query_directed_basic(matrix_species_names_a, matrix_a,
                                         matrix_species_names_b, matrix_b, 
                                         std::back_inserter(result_a_to_b),
                                         std::back_inserter(result_b_to_a));
      }   

    } // if(matrix_b.size() > 0)
    else
    {
      if(query_pairs.size() > 0)
      {
        cdnt.matrix_query_directed_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                                        query_pairs, std::back_inserter(result_a_to_b),
                                                        std::back_inserter(result_b_to_a)); 
      } 
      else
      {
        cdnt.matrix_query_directed_basic(matrix_species_names_a, matrix_a,
                                         std::back_inserter(result_a_to_b));
      }   

    } // else of if(matrix_b.size() > 0) 


    for(int i=0; i<result_a_to_b.size(); i++)
      output[i] = result_a_to_b[i];

    if(query_pairs.size() > 0 || matrix_b.size()>0 )
      for(int i=0; i<result_b_to_a.size(); i++)
        output[i+result_a_to_b.size()] = result_b_to_a[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void cdnt_query_directed(...)


void phylosor_query(int *n_w, int *n_l, int *edges, double *weights, char **names,
                    char **matrix_names_a, int *n_r_a, int *n_c_a, int *matrix_data_a, 
                    char **matrix_names_b, int *n_r_b, int *n_c_b, int *matrix_data_b,
                    int *n_p, int *pair_data, double *output,
                    char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names_a, matrix_species_names_b;
  std::vector<std::vector<bool> > matrix_a, matrix_b;
  std::vector<std::pair<int,int> > query_pairs;  

  transform_matrix_query_arguments_bimodal(n_w, n_l, edges, weights, names, 
                                           matrix_names_a, n_r_a, n_c_a, matrix_data_a,
                                           matrix_names_b, n_r_b, n_c_b, matrix_data_b,
                                           n_p, pair_data, parents, children, 
                                           edge_weights, species_names, 
                                           matrix_species_names_a, matrix_a,
                                           matrix_species_names_b, matrix_b, 
                                           query_pairs);

  try
  {
    Bimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Phylogenetic_Sorensens_similarity phylosor(tree);

    std::vector<Number_type> result;

    if(matrix_b.size() > 0)
    {
      if(query_pairs.size() > 0)
      {
        phylosor.matrix_query_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                                   matrix_species_names_b, matrix_b, 
                                                   query_pairs, std::back_inserter(result)); 
      }
      else
      {
        phylosor.matrix_query_basic(matrix_species_names_a, matrix_a,
                                    matrix_species_names_b, matrix_b, 
                                    std::back_inserter(result));
      }   

    } // if(matrix_b.size() > 0)
    else
    {
      if(query_pairs.size() > 0)
      {
        phylosor.matrix_query_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                                   query_pairs, std::back_inserter(result)); 
      } 
      else
      {
        phylosor.matrix_query_basic(matrix_species_names_a, matrix_a,
                                    std::back_inserter(result));
      }   

    } // else of if(matrix_b.size() > 0) 


    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void phylosor_query(...)


void unifrac_query(int *n_w, int *n_l, int *edges, double *weights, char **names,
                   char **matrix_names_a, int *n_r_a, int *n_c_a, int *matrix_data_a, 
                   char **matrix_names_b, int *n_r_b, int *n_c_b, int *matrix_data_b,
                   int *n_p, int *pair_data, double *output,
                   char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names_a, matrix_species_names_b;
  std::vector<std::vector<bool> > matrix_a, matrix_b;
  std::vector<std::pair<int,int> > query_pairs;  

  transform_matrix_query_arguments_bimodal(n_w, n_l, edges, weights, names, 
                                           matrix_names_a, n_r_a, n_c_a, matrix_data_a,
                                           matrix_names_b, n_r_b, n_c_b, matrix_data_b,
                                           n_p, pair_data, parents, children, 
                                           edge_weights, species_names, 
                                           matrix_species_names_a, matrix_a,
                                           matrix_species_names_b, matrix_b, 
                                           query_pairs);

  try
  {
    Bimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Unique_fraction unifrac(tree);

    std::vector<Number_type> result;

    if(matrix_b.size() > 0)
    {
      if(query_pairs.size() > 0)
      {
        unifrac.matrix_query_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                                  matrix_species_names_b, matrix_b, 
                                                  query_pairs, std::back_inserter(result)); 
      }
      else
      {
        unifrac.matrix_query_basic(matrix_species_names_a, matrix_a,
                                   matrix_species_names_b, matrix_b, 
                                   std::back_inserter(result));
      }   

    } // if(matrix_b.size() > 0)
    else
    {
      if(query_pairs.size() > 0)
      {
        unifrac.matrix_query_specific_pairs_basic(matrix_species_names_a, matrix_a,
                                                  query_pairs, std::back_inserter(result)); 
      } 
      else
      {
        unifrac.matrix_query_basic(matrix_species_names_a, matrix_a,
                                   std::back_inserter(result));
      }   

    } // else of if(matrix_b.size() > 0) 


    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void unifrac_query(...)

////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////////// Moments functions /////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

void pd_moments(int *n_w, int *n_l, int *edges, double *weights, char **names,
                int *n_q, int *q_sizes, bool *comp_expectation, bool *comp_deviation, 
                double *output, char **error_message, int *message_size )
{
  std::vector<int> parents, children, query_sizes;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names;

  transform_moments_function_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                                n_q, q_sizes, parents, children, 
                                                edge_weights, species_names, query_sizes);
  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Phylogenetic_diversity pd(tree);

    std::vector<Number_type> result;

    if(*comp_expectation == true)
      for(int i=0; i<query_sizes.size(); i++)
        output[i] = pd.compute_expectation(query_sizes[i]);  

    if(*comp_deviation == true)
    {
      if(*comp_expectation == true)
        for(int i=0; i<query_sizes.size(); i++)
          output[i+query_sizes.size()] = pd.compute_deviation(query_sizes[i]);
      else
        for(int i=0; i<query_sizes.size(); i++)
          output[i] = pd.compute_deviation(query_sizes[i]);
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void pd_moments(...)

void mpd_moments(int *n_w, int *n_l, int *edges, double *weights, char **names,
                 int *n_q, int *q_sizes, bool *comp_expectation, bool *comp_deviation, 
                 double *output, char **error_message, int *message_size )
{
  std::vector<int> parents, children, query_sizes;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names;

  transform_moments_function_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                                n_q, q_sizes, parents, children, 
                                                edge_weights, species_names, query_sizes);
  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Mean_pairwise_distance mpd(tree);

    std::vector<Number_type> result;

    if(*comp_expectation == true)
      for(int i=0; i<query_sizes.size(); i++)
        output[i] = mpd.compute_expectation(query_sizes[i]);  

    if(*comp_deviation == true)
    {
      if(*comp_expectation == true)
        for(int i=0; i<query_sizes.size(); i++)
          output[i+query_sizes.size()] = mpd.compute_deviation(query_sizes[i]); 
      else
        for(int i=0; i<query_sizes.size(); i++)
          output[i] = mpd.compute_deviation(query_sizes[i]);
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void mpd_moments(...)


void mntd_moments(int *n_w, int *n_l, int *edges, double *weights, char **names,
                  int *n_q, int *q_sizes, bool *comp_expectation, bool *comp_deviation, 
                  double *output, char **error_message, int *message_size )
{
  std::vector<int> parents, children, query_sizes;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names;

  transform_moments_function_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                                n_q, q_sizes, parents, children, 
                                                edge_weights, species_names, query_sizes);
  try
  {
    Mean_nearest_taxon_distance_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    if(!tree.is_ultrametric())
      return;

    Mean_nearest_taxon_distance mntd(tree);

    std::vector<Number_type> result;

    if(*comp_expectation == true)
      for(int i=0; i<query_sizes.size(); i++)
        output[i] = mntd.compute_expectation(query_sizes[i]);  

    if(*comp_deviation == true)
    {
      if(*comp_expectation == true)
        for(int i=0; i<query_sizes.size(); i++)
          output[i+query_sizes.size()] = mntd.compute_deviation(query_sizes[i]); 
      else
        for(int i=0; i<query_sizes.size(); i++)
          output[i] = mntd.compute_deviation(query_sizes[i]);
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");   
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();  
  }

  return;

} // void mntd_moments(...)


//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
////// Moments functions with abundance weights //////
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

void pd_moments_abundance_weighted(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                   int *n_q, int *q_sizes, 
                                   char **abundance_weights_nms, 
                                   double *abundance_weights_vals, 
                                   bool *comp_expectation, 
                                   bool *comp_deviation, double *output, 
                                   char **error_message, int *message_size )
{
  std::vector<int> parents, children, query_sizes;
  std::vector<Number_type> edge_weights, abundance_weights_values;
  std::vector<std::string> species_names, abundance_weights_names;

  transform_moments_function_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                                n_q, q_sizes, parents, children, 
                                                edge_weights, species_names, query_sizes);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Poisson_PD poisson_pd;

    std::vector<Number_type> exps, vars;

    int sample_size=0;

    for(int i=0; i<query_sizes.size(); i++)
      if(query_sizes[i] > sample_size)
        sample_size = query_sizes[i];

    bool level_based = true;

    poisson_pd.compute_expectations_and_variances(tree, sample_size, 
                                                  std::back_inserter(exps), std::back_inserter(vars), 
                                                  level_based);

    if(*comp_expectation == true)
      for(int i=0; i<query_sizes.size(); i++)
        output[i] = exps[query_sizes[i]];  

    if(*comp_deviation == true)
    {
      if(*comp_expectation == true)
        for(int i=0; i<query_sizes.size(); i++)
          output[i+query_sizes.size()] = Square_root()(vars[query_sizes[i]]);
      else
        for(int i=0; i<query_sizes.size(); i++)
          output[i] = Square_root()(vars[query_sizes[i]]);
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void pd_moments_abundance_weighted(...)


void mpd_moments_abundance_weighted(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                   int *n_q, int *q_sizes, 
                                   char **abundance_weights_nms, 
                                   double *abundance_weights_vals, 
                                   bool *comp_expectation, 
                                   bool *comp_deviation, double *output, 
                                   char **error_message, int *message_size )
{
  std::vector<int> parents, children, query_sizes;
  std::vector<Number_type> edge_weights, abundance_weights_values;
  std::vector<std::string> species_names, abundance_weights_names;

  transform_moments_function_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                                n_q, q_sizes, parents, children, 
                                                edge_weights, species_names, query_sizes);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Poisson_MPD poisson_mpd;

    std::vector<Number_type> exps, vars;

    int sample_size=0;

    for(int i=0; i<query_sizes.size(); i++)
      if(query_sizes[i] > sample_size)
        sample_size = query_sizes[i];

    bool level_based = true;

    poisson_mpd.compute_expectations_and_variances(tree, sample_size, 
                                                  std::back_inserter(exps), std::back_inserter(vars), 
                                                  level_based);

    if(*comp_expectation == true)
      for(int i=0; i<query_sizes.size(); i++)
        output[i] = exps[query_sizes[i]];  

    if(*comp_deviation == true)
    {
      if(*comp_expectation == true)
        for(int i=0; i<query_sizes.size(); i++)
          output[i+query_sizes.size()] = Square_root()(vars[query_sizes[i]]);
      else
        for(int i=0; i<query_sizes.size(); i++)
          output[i] = Square_root()(vars[query_sizes[i]]);
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void mpd_moments_abundance_weighted(...)


void mntd_moments_abundance_weighted(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                     int *n_q, int *q_sizes, 
                                     char **abundance_weights_nms, 
                                     double *abundance_weights_vals, 
                                     bool *comp_expectation, 
                                     bool *comp_deviation, double *output, 
                                     char **error_message, int *message_size )
{
  std::vector<int> parents, children, query_sizes;
  std::vector<Number_type> edge_weights, abundance_weights_values;
  std::vector<std::string> species_names, abundance_weights_names;

  transform_moments_function_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                                n_q, q_sizes, parents, children, 
                                                edge_weights, species_names, query_sizes);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Mean_nearest_taxon_distance_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Poisson_MNTD poisson_mntd;

    std::vector<Number_type> exps, vars;

    int sample_size=0;

    for(int i=0; i<query_sizes.size(); i++)
      if(query_sizes[i] > sample_size)
        sample_size = query_sizes[i];

    bool level_based = true;

    poisson_mntd.compute_expectations_and_variances(tree, sample_size, 
                                                    std::back_inserter(exps), std::back_inserter(vars), 
                                                    level_based);

    if(*comp_expectation == true)
      for(int i=0; i<query_sizes.size(); i++)
        output[i] = exps[query_sizes[i]];  

    if(*comp_deviation == true)
    {
      if(*comp_expectation == true)
        for(int i=0; i<query_sizes.size(); i++)
          output[i+query_sizes.size()] = Square_root()(vars[query_sizes[i]]);
      else
        for(int i=0; i<query_sizes.size(); i++)
          output[i] = Square_root()(vars[query_sizes[i]]);
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void mntd_moments_abundance_weighted(...)

///////////////////////////////////////////////////////////////////////////////
// Computation of weighted moments using the R-style sequential distribution //
///////////////////////////////////////////////////////////////////////////////

void pd_moments_weighted_sequential(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                    int *n_q, int *q_sizes, 
                                    char **abundance_weights_nms, 
                                    double *abundance_weights_vals, 
                                    bool *comp_expectation, 
                                    bool *comp_deviation, int *repetitions, int *seed, 
                                    double *output, char **error_message, int *message_size )
{
  std::vector<int> parents, children, query_sizes;
  std::vector<Number_type> edge_weights, abundance_weights_values;
  std::vector<std::string> species_names, abundance_weights_names;

  transform_moments_function_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                                n_q, q_sizes, parents, children, 
                                                edge_weights, species_names, query_sizes);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                              abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);
 
    std::vector< std::pair<Number_type,Number_type> > moments;

    Phylogenetic_diversity pd(tree);

    std::vector<int> leaf_indices;
    std::vector<Number_type> probs;
    typename Unimodal_tree::Leaves_iterator it;
 
    for(it = tree.leaves_begin(); it != tree.leaves_end(); it++)
    {
      leaf_indices.push_back(it->second);
      probs.push_back(tree.node_probability(it->second));
    }  
    
    Sequential_sampler weighted_sampler(leaf_indices,probs);
    
    pd.set_seed(*seed);

    Incremental_Monte_Carlo_handler MC_handler;

    MC_handler.estimate_moments_with_Monte_Carlo(pd, query_sizes, weighted_sampler, *repetitions, std::back_inserter(moments));


    if(*comp_expectation == true)
      for(int i=0; i<moments.size(); i++)
        output[i] = moments[i].first;  

    if(*comp_deviation == true)
    {
      if(*comp_expectation == true)
        for(int i=0; i<query_sizes.size(); i++)
          output[i+query_sizes.size()] = moments[i].second;
      else
        for(int i=0; i<query_sizes.size(); i++)
          output[i] = moments[i].second;
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void pd_moments_weighted_sequential(...)


void mpd_moments_weighted_sequential(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                     int *n_q, int *q_sizes, 
                                     char **abundance_weights_nms, 
                                     double *abundance_weights_vals, 
                                     bool *comp_expectation, 
                                     bool *comp_deviation, int *repetitions, int *seed, 
                                     double *output, char **error_message, int *message_size )
{
  std::vector<int> parents, children, query_sizes;
  std::vector<Number_type> edge_weights, abundance_weights_values;
  std::vector<std::string> species_names, abundance_weights_names;

  transform_moments_function_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                                n_q, q_sizes, parents, children, 
                                                edge_weights, species_names, query_sizes);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                              abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);
 
    std::vector< std::pair<Number_type,Number_type> > moments;

    Mean_pairwise_distance mpd(tree);

    std::vector<int> leaf_indices;
    std::vector<Number_type> probs;
    typename Unimodal_tree::Leaves_iterator it;
 
    for(it = tree.leaves_begin(); it != tree.leaves_end(); it++)
    {
      leaf_indices.push_back(it->second);
      probs.push_back(tree.node_probability(it->second));
    }  
    
    Sequential_sampler weighted_sampler(leaf_indices,probs);

    mpd.set_seed(*seed);

    Incremental_Monte_Carlo_handler MC_handler;

    MC_handler.estimate_moments_with_Monte_Carlo(mpd, query_sizes, weighted_sampler, *repetitions, std::back_inserter(moments));


    if(*comp_expectation == true)
      for(int i=0; i<moments.size(); i++)
        output[i] = moments[i].first;  

    if(*comp_deviation == true)
    {
      if(*comp_expectation == true)
        for(int i=0; i<query_sizes.size(); i++)
          output[i+query_sizes.size()] = moments[i].second;
      else
        for(int i=0; i<query_sizes.size(); i++)
          output[i] = moments[i].second;
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void mpd_moments_weighted_sequential(...)


void mntd_moments_weighted_sequential(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                      int *n_q, int *q_sizes, 
                                      char **abundance_weights_nms, 
                                      double *abundance_weights_vals, 
                                      bool *comp_expectation, 
                                      bool *comp_deviation, int *repetitions, int *seed, 
                                      double *output, char **error_message, int *message_size )
{
  std::vector<int> parents, children, query_sizes;
  std::vector<Number_type> edge_weights, abundance_weights_values;
  std::vector<std::string> species_names, abundance_weights_names;

  transform_moments_function_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                                n_q, q_sizes, parents, children, 
                                                edge_weights, species_names, query_sizes);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                              abundance_weights_names, abundance_weights_values);

  try
  {
    Mean_nearest_taxon_distance_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);
 
    std::vector< std::pair<Number_type,Number_type> > moments;

    Mean_nearest_taxon_distance mntd(tree);

    std::vector<int> leaf_indices;
    std::vector<Number_type> probs;
    typename Mean_nearest_taxon_distance_tree::Leaves_iterator it;
 
    for(it = tree.leaves_begin(); it != tree.leaves_end(); it++)
    {
      leaf_indices.push_back(it->second);
      probs.push_back(tree.node_probability(it->second));
    }  
    
    Sequential_sampler weighted_sampler(leaf_indices,probs);

    mntd.set_seed(*seed);

    Incremental_Monte_Carlo_handler MC_handler;

    MC_handler.estimate_moments_with_Monte_Carlo(mntd, query_sizes, weighted_sampler, *repetitions, std::back_inserter(moments));


    if(*comp_expectation == true)
      for(int i=0; i<moments.size(); i++)
        output[i] = moments[i].first;  

    if(*comp_deviation == true)
    {
      if(*comp_expectation == true)
        for(int i=0; i<query_sizes.size(); i++)
          output[i+query_sizes.size()] = moments[i].second;
      else
        for(int i=0; i<query_sizes.size(); i++)
          output[i] = moments[i].second;
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void mntd_moments_weighted_sequential(...)

void cac_moments_weighted_sequential(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                     double *chi, int *n_q, int *q_sizes, 
                                     char **abundance_weights_nms, 
                                     double *abundance_weights_vals, 
                                     bool *comp_expectation, 
                                     bool *comp_deviation, int *repetitions, int *seed, 
                                     double *output, char **error_message, int *message_size )
{
  std::vector<int> parents, children, query_sizes;
  std::vector<Number_type> edge_weights, abundance_weights_values;
  std::vector<std::string> species_names, abundance_weights_names;

  transform_moments_function_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                                n_q, q_sizes, parents, children, 
                                                edge_weights, species_names, query_sizes);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                              abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);
 
    std::vector< std::pair<Number_type,Number_type> > moments;

    Core_ancestor_cost cac(tree,*chi);

    std::vector<int> leaf_indices;
    std::vector<Number_type> probs;
    typename Unimodal_tree::Leaves_iterator it;
 
    for(it = tree.leaves_begin(); it != tree.leaves_end(); it++)
    {
      leaf_indices.push_back(it->second);
      probs.push_back(tree.node_probability(it->second));
    }  
    
    Sequential_sampler weighted_sampler(leaf_indices,probs);

    cac.set_seed(*seed);

    Incremental_Monte_Carlo_handler MC_handler;

    MC_handler.estimate_moments_with_Monte_Carlo(cac, query_sizes, weighted_sampler, *repetitions, std::back_inserter(moments));


    if(*comp_expectation == true)
      for(int i=0; i<moments.size(); i++)
        output[i] = moments[i].first;  

    if(*comp_deviation == true)
    {
      if(*comp_expectation == true)
        for(int i=0; i<query_sizes.size(); i++)
          output[i+query_sizes.size()] = moments[i].second;
      else
        for(int i=0; i<query_sizes.size(); i++)
          output[i] = moments[i].second;
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void cac_moments_weighted_sequential(...)

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

void cac_moments(int *n_w, int *n_l, int *edges, double *weights, char **names, double *chi,
                 int *n_q, int *q_sizes, int *k, double *output,
                 char **error_message, int *message_size )
{
  std::vector<int> parents, children, query_sizes;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names;

  transform_moments_function_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                                n_q, q_sizes, parents, children, 
                                                edge_weights, species_names, query_sizes);
  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    if(!tree.is_ultrametric())
      return;

    Core_ancestor_cost cac(tree, *chi);

    for(int i=0; i<query_sizes.size(); i++)
    {
      std::vector<Number_type> result;
      cac.compute_first_k_centralised_moments( *k, query_sizes[i], std::back_inserter(result));

      for(int j=0; j<*k; j++)
        output[(j*query_sizes.size()) + i] = result[j]; 
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");   
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void cac_moments(...)


void cbl_moments(int *n_w, int *n_l, int *edges, double *weights, char **names,
                 int *n_q, int *q_sizes, bool *comp_expectation, bool *comp_deviation, 
                 double *output, char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names;
  std::vector< std::pair<int,int> > query_sizes;

  transform_moments_function_arguments_bimodal(n_w, n_l, edges, weights, names, 
                                               n_q, q_sizes, parents, children, 
                                               edge_weights, species_names, query_sizes);
  try
  {
    Bimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Common_branch_length cbl(tree);

    std::vector<Number_type> result;

    if(*comp_expectation == true)
      for(int i=0; i<query_sizes.size(); i++)
        output[i] = cbl.compute_expectation(query_sizes[i].first, query_sizes[i].second);  

    if(*comp_deviation == true)
    {
      if(*comp_expectation == true)
        for(int i=0; i<query_sizes.size(); i++)
          output[i+query_sizes.size()] = cbl.compute_deviation(query_sizes[i].first, query_sizes[i].second); 
      else
        for(int i=0; i<query_sizes.size(); i++)
          output[i] = cbl.compute_deviation(query_sizes[i].first, query_sizes[i].second);
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");   
    copy_message(error_msg,error_message);
    *message_size = error_msg.size(); 
  }

  return;

} // void cbl_moments(...)

void cd_moments(int *n_w, int *n_l, int *edges, double *weights, char **names,
                int *n_q, int *q_sizes, bool *comp_expectation, bool *comp_deviation, 
                double *output, char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names;
  std::vector< std::pair<int,int> > query_sizes;

  transform_moments_function_arguments_bimodal(n_w, n_l, edges, weights, names, 
                                               n_q, q_sizes, parents, children, 
                                               edge_weights, species_names, query_sizes);
  try
  {
    Bimodal_tree tree;
  
    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    Community_distance cd(tree);

    std::vector<Number_type> result;

    if(*comp_expectation == true)
      for(int i=0; i<query_sizes.size(); i++)
        output[i] = cd.compute_expectation(query_sizes[i].first, query_sizes[i].second);  

    if(*comp_deviation == true)
    {
      if(*comp_expectation == true)
        for(int i=0; i<query_sizes.size(); i++)
          output[i+query_sizes.size()] = cd.compute_deviation(query_sizes[i].first, query_sizes[i].second); 
      else
        for(int i=0; i<query_sizes.size(); i++)
          output[i] = cd.compute_deviation(query_sizes[i].first, query_sizes[i].second);
    }

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");   
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();  
  }

  return;

} // void cd_moments(...)

///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////
/// Functions that compute p-values ///
///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////// Uniform p-value functions that use sequential Monte-Carlo computations ///////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void pd_pvalues_uniform(int *n_w, int *n_l, int *edges, double *weights, char **names,
                         char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                         int *repetitions, int *seed, double *output,
                         char **error_message, int *message_size)
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    std::vector<Number_type> result;
    int reps = *repetitions; 

    Phylogenetic_diversity pd(tree);

    pd.set_probability_distribution(Kernel::UNIFORM_FIXED_SIZE);
    pd.set_seed(*seed);

    pd.pvalues_query_uniform_fixed_size(matrix_species_names, matrix, 
                                        std::back_inserter(result), reps);  
   
    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void pd_pvalues_uniform(...)


void mpd_pvalues_uniform(int *n_w, int *n_l, int *edges, double *weights, char **names,
                         char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                         int *repetitions, int *seed, double *output,
                         char **error_message, int *message_size)
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    std::vector<Number_type> result;
    int reps = *repetitions; 

    Mean_pairwise_distance mpd(tree);

    mpd.set_probability_distribution(Kernel::UNIFORM_FIXED_SIZE);

    mpd.set_seed(*seed);
    mpd.pvalues_query_uniform_fixed_size(matrix_species_names, 
                                         matrix, std::back_inserter(result),reps);    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void mpd_pvalues_uniform(...)

void mntd_pvalues_uniform(int *n_w, int *n_l, int *edges, double *weights, char **names,
                           char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                           int *repetitions, int *seed, double *output,
                           char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  try
  {
    Mean_nearest_taxon_distance_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    std::vector<Number_type> result;
    int reps = *repetitions;

    Mean_nearest_taxon_distance mntd(tree);

    mntd.set_probability_distribution(Kernel::UNIFORM_FIXED_SIZE);

    mntd.set_seed(*seed);
    mntd.pvalues_query_uniform_fixed_size(matrix_species_names, matrix, 
                                          std::back_inserter(result),reps);  
  
    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void mntd_pvalues_uniform(...)


void cac_pvalues_uniform(int *n_w, int *n_l, int *edges, double *weights,  
                         char **names, double *chi, char **matrix_names, 
                         int *n_r, int *n_c, int *matrix_data, 
                         int *repetitions, int *seed, double *output,
                         char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights;
  std::vector<std::string> species_names, matrix_species_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);

    std::vector<Number_type> result;
    int reps = *repetitions;

    Core_ancestor_cost cac(tree,*chi);

    cac.set_probability_distribution(Kernel::UNIFORM_FIXED_SIZE);
    cac.set_seed(*seed);

    cac.pvalues_query_uniform_fixed_size(matrix_species_names, matrix, 
                                         std::back_inserter(result),reps);  
    
    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void cac_pvalues_uniform(...)




/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////// Abundance-weighted p-value functions that use sequential Monte-Carlo computations ///////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void pd_pvalues_weighted_sequential(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                    char **abundance_weights_nms, double *abundance_weights_vals, 
                                    char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                                    int *repetitions, int *seed, double *output,
                                    char **error_message, int *message_size)
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights,abundance_weights_values;
  std::vector<std::string> species_names, matrix_species_names, abundance_weights_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Phylogenetic_diversity pd(tree);

    std::vector<Number_type> result;
    int reps = *repetitions; 

    pd.set_probability_distribution(Kernel::SEQUENTIAL_FIXED_SIZE);
    pd.set_seed(*seed);

    pd.pvalues_query_sequential_fixed_size(matrix_species_names, matrix, 
                                           std::back_inserter(result), reps);  
   

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void pd_pvalues_weighted_sequential(...)


void mpd_pvalues_weighted_sequential(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                     char **abundance_weights_nms, double *abundance_weights_vals, 
                                     char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                                     int *repetitions, int *seed, double *output,
                                     char **error_message, int *message_size)
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights,abundance_weights_values;
  std::vector<std::string> species_names, matrix_species_names, abundance_weights_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Mean_pairwise_distance mpd(tree);

    std::vector<Number_type> result;
    int reps = *repetitions; 

    mpd.set_probability_distribution(Kernel::SEQUENTIAL_FIXED_SIZE);
    mpd.set_seed(*seed);

    mpd.pvalues_query_sequential_fixed_size(matrix_species_names, 
                                            matrix, std::back_inserter(result),reps);    

    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message();
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void mpd_pvalues_weighted_sequential(...)

void mntd_pvalues_weighted_sequential(int *n_w, int *n_l, int *edges, double *weights, char **names,
                                      char **abundance_weights_nms, double *abundance_weights_vals, 
                                      char **matrix_names, int *n_r, int *n_c, int *matrix_data, 
                                      int *repetitions, int *seed, double *output,
                                      char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights, abundance_weights_values;
  std::vector<std::string> species_names, matrix_species_names,abundance_weights_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Mean_nearest_taxon_distance_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Mean_nearest_taxon_distance mntd(tree);

    mntd.set_probability_distribution(Kernel::SEQUENTIAL_FIXED_SIZE);

    std::vector<Number_type> result;
    int reps = *repetitions;

     mntd.set_seed(*seed);
     mntd.pvalues_query_sequential_fixed_size(matrix_species_names, matrix, 
                                              std::back_inserter(result),reps);  
  
    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void mntd_pvalues_weighted_sequential(...)


void cac_pvalues_weighted_sequential(int *n_w, int *n_l, int *edges, double *weights,  
                                     char **names, double *chi, char **abundance_weights_nms, 
                                     double *abundance_weights_vals, char **matrix_names, 
                                     int *n_r, int *n_c, int *matrix_data, 
                                     int *repetitions, int *seed, double *output,
                                     char **error_message, int *message_size )
{
  std::vector<int> parents, children;
  std::vector<Number_type> edge_weights, abundance_weights_values;
  std::vector<std::string> species_names, matrix_species_names,abundance_weights_names;
  std::vector<std::vector<bool> > matrix;

  transform_matrix_query_arguments_unimodal(n_w, n_l, edges, weights, names, 
                                            matrix_names, n_r, n_c, matrix_data,
                                            parents, children, edge_weights, 
                                            species_names, matrix_species_names, matrix);

  transform_abundance_weights(n_l, abundance_weights_nms, abundance_weights_vals,
                            abundance_weights_names, abundance_weights_values);

  try
  {
    Unimodal_tree tree;

    tree.construct_from_edge_data(parents, children, edge_weights, species_names);
    tree.set_leaf_probability_values(abundance_weights_names, abundance_weights_values);

    Core_ancestor_cost cac(tree,*chi);

    cac.set_probability_distribution(Kernel::SEQUENTIAL_FIXED_SIZE);

    std::vector<Number_type> result;
    int reps = *repetitions;

    cac.set_seed(*seed);
    cac.pvalues_query_sequential_fixed_size(matrix_species_names, matrix, 
                                            std::back_inserter(result),reps);  
    
    for(int i=0; i<result.size(); i++)
      output[i] = result[i];

    tree.clear();
    flush_warnings();
    *message_size = 0;

  } // end of try
  catch(Exception_type expc)
  {
    flush_warnings();

    std::string error_msg = expc.return_error_message(); 
    REprintf("\n");  
    copy_message(error_msg,error_message);
    *message_size = error_msg.size();
  }

  return;

} // void cac_pvalues_weighted_sequential(...)


} // extern "C"
