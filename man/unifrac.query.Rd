\name{unifrac.query}
\alias{unifrac.query}
\title{
Computes the value of the Unique Fraction measure
}
\description{
Calculates the Unique Fraction (UniFrac) given paired sets of tips on a phylogeny. 
A version of this function that computes  the standardised value of the measure is not yet available.
}
\usage{
unifrac.query(tree, matrix.a, matrix.b = NULL, 
              query.matrix = NULL)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{matrix.a}{A matrix with binary (0/1) values, where each row represents a tip set. Each column name in the matrix 
                  must match a tip label on the input tree. If not all values in the matrix are binary, 
                  we consider two cases; if the matrix contains only non-negative values, all values are 
                  coerced to binary ones and a warning message is printed. If the matrix contains at least one 
                  negative value, the function throws an error}
  \item{matrix.b}{Optional, a second matrix with a similar format as matrix.a}
  \item{query.matrix}{Optional, a two-column matrix specifying the pairs of rows (tip sets)
                      for which the function computes the UniFrac values. Each row in query.matrix
                      indicates a pair of tip sets for which we want to compute the UniFrac value.
                      Let k and r be the values that are stored in the i-th row of query.matrix, where k
                      is the value stored in the first column and r is the value stored 
                      in the second column. If matrix.b is given, the function computes the UniFrac value between the 
                      k-th row of matrix.a and the r-th row of matrix.b. If matrix.b is not given,
                      the function computes the UniFrac value between the k-th and r-th row of matrix.a (default = NULL)}
}
\details{Queries can be given in four ways. If neither matrix.b nor query.matrix are given, the function computes the UniFrac values for all pairs of rows (tip sets) in matrix.a . If matrix.b is given but not query.matrix, the function computes the UniFrac values for all combinations of a row in matrix.a with rows in matrix.b. If query.matrix is given and matrix.b is not, the function returns the UniFrac values for the pairs of rows in matrix.a specified by query.matrix. If query.matrix and matrix.b are both given, UniFrac values are computed for the rows in matrix.a specified by the first column of query.matrix against the rows in matrix.b specified in the second column of query.matrix.}
\value{The UniFrac values for the requested pairs of tip sets. If query.matrix is provided, then the values are returned in an one-dimensional vector. The i-th element of this vector is the UniFrac value for the pair of tip sets indicated in the i-th row of query.matrix. If query.matrix is not provided, the UniFrac values are returned in a matrix object; entry [i,j] in the output matrix stores the UniFrac value between the tip sets specified on the i-th and j-th row of matrix.a (if matrix.b is not specified), or the UniFrac value between the i-th row of matrix.a and the j-th row of matrix.b (if matrix.b is specified) }
\references{
Lozupone C. and R. Knight. 2005. UniFrac: a New Phylogenetic Method for Comparing Microbial Communities. Applied and Environmental Microbiology,71(12):8228-35.}

\author{Constantinos Tsirogiannis (tsirogiannis.c@gmail.com) and Brody Sandel (bsandel@scu.edu)}

\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

#Create 10 random communities with 50 families each
comm = matrix(0,nrow = 10,ncol = length(bird.families$tip.label))
for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
colnames(comm) = bird.families$tip.label

#Calculate all pairwise UniFrac values for communities in comm
unifrac.query(bird.families,comm)

#Calculate pairwise distances from 
#the first two rows of comm to all rows
unifrac.query(bird.families, comm[1:2,],comm)

#Calculate the distances from the first two rows 
#to all rows using the query matrix
qm = expand.grid(1:2,1:10)
unifrac.query(bird.families,comm,query.matrix = qm)     
}

