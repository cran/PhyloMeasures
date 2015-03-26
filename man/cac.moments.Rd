\name{cac.moments}
\alias{cac.moments}
\title{
Computes the statistical moments of the Core Ancestor Cost measure
}
\description{
Calculates the first k statistical moments of the Core Ancestor Cost measure for a tree, a vector of tip set sizes, and a value of chi. Contrary to other measures, this function makes it possible to compute moments which are of higher order than the mean and the variance. The CAC is calculated by identifying the node in the tree that is the most recent common ancestor of at least chi proportion of the tips in the set. The CAC is the distance of this node from the root of the tree.
}
\usage{
cac.moments(tree, chi, sample.sizes, k=2)
}
\arguments{
  \item{tree}{A phylo tree object}
  \item{chi}{A number in the interval (0.5,1]}
  \item{sample.sizes}{A vector of non-negative integers specifying the tip set sizes for which to calculate moments}
  \item{k}{A positive integer specifying the number of moments to compute (default = 2)}

}
\value{A k-column matrix with length(sample.sizes) rows. Entry [i,j]  in the matrix gives the j-th moment for the i-th sample size in sample.sizes . The first moment in each row is the mean, and for j larger than one, the j-th returned moment is the centralised statistical moment of order j. Given a phylogenetic tree T, a sample size r and a value chi, the 
centralised statistical moment of order j for the CAC is defined as: 
\deqn{E_{R \in Sub(T,r)}[(\ensuremath{\mathrm{CAC}}_{\chi}(T,R)-m_{\ensuremath{\mathrm{CAC}}}(T,r))^k],}{E_{R in Sub(T,r)}[(CAC_{chi}(T,R)-m_{CAC,chi}(T,r))^k],}
where \eqn{E_{R \in Sub(T,r)}}{E_{R \in Sub(T,r)}} denotes the expectation of a random variable for all tip sets
in T that consist of r tips each, and \eqn{m_{\ensuremath{\mathrm{CAC}}}(T,r)}{m_{CAC,chi}(T,r)} 
indicates the mean value of the CAC over all these sets.}
\references{Tsirogiannis, C., B. Sandel and A. Kalvisa. 2014. New algorithms for computing phylogenetic biodiversity. Algorithms in Bioinformatics, LNCS 8701: 187-203.}
\author{Constantinos Tsirogiannis (constant@madalgo.au.dk)}

\seealso{\code{\link{cac.query}}
}
\examples{
#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

cac.moments(bird.families,0.75,1:100,k=4)
}
