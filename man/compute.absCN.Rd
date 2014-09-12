\name{compute.absCN}
\alias{compute.absCN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~ function to do. ~~
Convert copy ratios into abosolute copy numbers when given the purity and ploidy estimates
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function converts copy ratios into abosolute copy numbers when the tumor purity and ploidy are provided by the user. 
The tumor purity and ploidy parameter should be taken from the result after running the "run.absCNSeq" function. 
}
\usage{
compute.absCN(seg.data, alpha.opt, tau.opt, qmax = 7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg.data}{
%%     ~~Describe \code{seg.data} here~~
a data frame with five columns: "chrom", "loc.start", "loc.end", "eff.seg.len", "normalized.ratio".
}
  \item{alpha.opt}{
%%     ~~Describe \code{alpha.opt} here~~
tumor purity estimate
}
  \item{tau.opt}{
%%     ~~Describe \code{tau.opt} here~~
tumor ploidy estimate
}
  \item{qmax}{
%%     ~~Describe \code{qmax} here~~
maximum allowed absolute copy number for any segments
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The input data frame is augmented with two additional columns: rhat (expected copy ratio) and CN (absolute copy number)
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Lei Bao
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
%% seg.data <- read.table("example.cn.txt",header=TRUE,sep="\t")
%% data(absCNseq)
%% seg.CN <- compute.absCN(seg.data, 0.58, 2.18)

library(absCNseq)
my.res.list <- run.absCNSeq("example.cn.txt", "example.snv.txt", "myResult", "Sample1", seq.type="WES", min.seg.len=200)
seg.CN <- compute.absCN(my.res.list$seg.dat, my.res.list$searchRes[i,"alpha"], my.res.list$searchRes[i,"tau"])  # the i-th solution

% ## The function is currently defined as
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
