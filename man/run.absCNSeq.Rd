\name{run.absCNSeq}
\alias{run.absCNSeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~ function to do. ~~
main function to run absCNseq 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This wrapper function accepts data files and user specified parameters and runs the absCNseq algorithm. 
}
\usage{
run.absCNSeq(seg.fn, snv.fn = NULL, res.dir, smp.name, seq.type = c("WES", "WGS"), alpha.min = 0.2, alpha.max = 1, tau.min = 1.5, tau.max = 5, min.freq = 0, min.seg.len = 0, qmax = 7, lamda = 0.5, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg.fn}{
%%     ~~Describe \code{seg.fn} here~~
The name of the file containing the segmentation data. An example file called "example.cn.txt" can be found under the "example" folder within the tarball.   
}
  \item{snv.fn}{
%%     ~~Describe \code{snv.fn} here~~
The name of the file containing a set of somatic single nucleotide variants (SNVs). An example file called "example.snv.txt" can be found under the "example" folder within the tarball.   
}
  \item{res.dir}{
%%     ~~Describe \code{res.dir} here~~
The output directory 
}
  \item{smp.name}{
%%     ~~Describe \code{smp.name} here~~
Sample name
}
  \item{seq.type}{
%%     ~~Describe \code{seq.type} here~~
Either "WES" (whole exome sequencing) or "WGS" (whole genome sequencing)
}
  \item{alpha.min}{
%%     ~~Describe \code{alpha.min} here~~
The minimum allowed value for tumor purity. Default is 0.20. If you do have the pathologist estimate, set it as the lower bound of the pathologist estimate is usually preferred. 
}
  \item{alpha.max}{
%%     ~~Describe \code{alpha.max} here~~
The maximum allowed value for tumor purity. Default is 1.0. If you do have the pathologist estimate, set it as the upper bound of the pathologist estimate is usually preferred. 
}
  \item{tau.min}{
%%     ~~Describe \code{tau.min} here~~
The minimum allowed value for tumor ploidy
}
  \item{tau.max}{
%%     ~~Describe \code{tau.max} here~~
The maximum allowed value for tumor ploidy
}
  \item{min.sol.freq}{
%%     ~~Describe \code{min.sol.freq} here~~
A solution should appear at least this many times to be kept. Singleton solutions are usually not trustable. By default (min.sol.freq=0), the program will only retain solutions that cover at least 1 percent of the search space.
}
  \item{min.seg.len}{
%%     ~~Describe \code{min.seg.len} here~~
The minimum length of a segment to be included in computation. The default value is 200 bp for WES  and 3000 bp for WGS. 
}
  \item{qmax}{
%%     ~~Describe \code{qmax} here~~
Maximum allowed absolute copy number for any segments.
}
  \item{lamda}{
%%     ~~Describe \code{lamda} here~~
The relative weight of the segment copy ratio data over the SNV data. Must be a value in(0.0,1.0]. 
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Please refer to the "example" folder to see the format of the input segmentation file or SNV file.
  
}
\value{
%%  ~Describe the value returned
a list is returned
%%  If it is a LIST, use
  \item{searchRes }{a data frame giving the solution set (purity and ploidy pairs) ranked by their fitting errors }
  \item{absCN }{a data frame giving the absolute copy number estimates of tumor cells for the top first solution (purity and ploidy pair)}
  \item{absSNV }{a data frame giving the absolute multiplicity estimates of SNVs for the top first solution (purity and ploidy pair) }
  \item{orig.seg.dat }{original copy ratio data read from the input segmentation file } 
  \item{seg.dat }{filtered copy ratio data } 
  \item{orig.snv.dat }{original SNV data read from the input SNV file }  
  \item{snv.dat }{filtered SNV data }
  
%% ...
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

library(absCNseq)
my.res.list <- run.absCNSeq("example.cn.txt", "example.snv.txt", "myResult", "Sample1", seq.type="WES", min.seg.len=200)
seg.CN <- compute.absCN(my.res.list$seg.dat, my.res.list$searchRes[i,"alpha"], my.res.list$searchRes[i,"tau"])  # the i-th solution
plot.absCN(seg.CN, chromnum=1)  # plot chromosome 1

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ main function }
%\keyword{ absCNSeq } % __ONLY ONE__ keyword per line
}