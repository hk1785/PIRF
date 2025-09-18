\name{immunotherapy}
\alias{immunotherapy}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("immunotherapy")}
\format{
  The format is:
Formal class 'phyloseq' [package "phyloseq"] with 5 slots
  ..@ otu_table:Formal class 'otu_table' [package "phyloseq"] with 2 slots
  .. .. ..@ .Data        : num [1:257, 1:986] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:257] "S1" "S2" "S3" "S4" ...
  .. .. .. .. ..$ : chr [1:986] "A1" "A2" "A3" "A4" ...
  .. .. ..@ taxa_are_rows: logi FALSE
  .. .. ..$ dim     : int [1:2] 257 986
  .. .. ..$ dimnames:List of 2
  .. .. .. ..$ : chr [1:257] "S1" "S2" "S3" "S4" ...
  .. .. .. ..$ : chr [1:986] "A1" "A2" "A3" "A4" ...
  ..@ tax_table: NULL
  ..@ sam_data :'data.frame':	257 obs. of  1 variable:
Formal class 'sample_data' [package "phyloseq"] with 4 slots
  .. .. ..@ .Data    :List of 1
  .. .. .. ..$ : num [1:257] 0 1 0 1 0 1 1 0 0 0 ...
  .. .. ..@ names    : chr "y"
  .. .. ..@ row.names: chr [1:257] "S1" "S2" "S3" "S4" ...
  .. .. ..@ .S3Class : chr "data.frame"
  ..@ phy_tree :List of 4
  .. ..$ edge       : int [1:1970, 1:2] 987 988 989 990 991 992 993 994 995 996 ...
  .. ..$ edge.length: num [1:1970] 0.0997 0.0749 0.3624 0.3568 0.4026 ...
  .. ..$ Nnode      : int 985
  .. ..$ tip.label  : chr [1:986] "A1" "A2" "A3" "A4" ...
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
data(immunotherapy)
## maybe str(immunotherapy) ; plot(immunotherapy) ...
}
\keyword{datasets}
