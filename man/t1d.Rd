\name{t1d}
\alias{t1d}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("t1d")}
\format{
  The format is:
Formal class 'phyloseq' [package "phyloseq"] with 5 slots
  ..@ otu_table:Formal class 'otu_table' [package "phyloseq"] with 2 slots
  .. .. ..@ .Data        : num [1:521, 1:307] 0.001443 0.000282 0 0 0.000294 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:521] "S1" "S2" "S3" "S4" ...
  .. .. .. .. ..$ : chr [1:307] "A1" "A2" "A3" "A4" ...
  .. .. ..@ taxa_are_rows: logi FALSE
  .. .. ..$ dim     : int [1:2] 521 307
  .. .. ..$ dimnames:List of 2
  .. .. .. ..$ : chr [1:521] "S1" "S2" "S3" "S4" ...
  .. .. .. ..$ : chr [1:307] "A1" "A2" "A3" "A4" ...
  ..@ tax_table: NULL
  ..@ sam_data :'data.frame':	521 obs. of  1 variable:
Formal class 'sample_data' [package "phyloseq"] with 4 slots
  .. .. ..@ .Data    :List of 1
  .. .. .. ..$ : int [1:521] 1 0 1 1 0 1 1 1 0 0 ...
  .. .. ..@ names    : chr "y"
  .. .. ..@ row.names: chr [1:521] "S1" "S2" "S3" "S4" ...
  .. .. ..@ .S3Class : chr "data.frame"
  ..@ phy_tree :List of 5
  .. ..$ edge       : int [1:612, 1:2] 308 309 310 311 312 313 314 315 315 316 ...
  .. ..$ edge.length: num [1:612] 0.00345 0.14015 0.08422 0.025 0.08719 ...
  .. ..$ Nnode      : int 306
  .. ..$ node.label : chr [1:306] "" "0.910" "1.000" "0.999" ...
  .. ..$ tip.label  : chr [1:307] "A1" "A2" "A3" "A4" ...
  .. ..- attr(*, "class")= chr "phylo"
  .. ..- attr(*, "order")= chr "cladewise"
  ..@ refseq   : NULL
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(t1d)
## maybe str(t1d) ; plot(t1d) ...
}
\keyword{datasets}
