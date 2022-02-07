\name{sample_data}
\alias{sample_data}
\docType{data}
\title{sample data for analysis of  NBB/PNBB design
}
\description{
Data is given here for illustration purpose.Here the data frame contains 5 columns of block, treatment, left_neighbour, right_neighbour and yield.
}
\usage{
sample_data
}
\format{
  A data frame with 110 observations on the following 5 variables.
  \describe{
    \item{\code{block}}{a factor with levels \code{1}  \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10}}
    \item{\code{treatment}}{a factor with levels \code{1}  \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10} \code{11}}
    \item{\code{left_neighbour}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10} \code{11}}
    \item{\code{right_neighbour}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10} \code{11}}
    \item{\code{yield}}{a numeric vector}
  }
}

\examples{
#To view the sample data set
library(NBBDesigns)
sample_data
#To analyze the sample_data
anlys(sample_data)
}

