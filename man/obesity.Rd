\name{obesity}
\alias{obesity}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("obesity")}
\format{
  The format is:
Formal class 'phyloseq' [package "phyloseq"] with 5 slots
  ..@ otu_table:Formal class 'otu_table' [package "phyloseq"] with 2 slots
  .. .. ..@ .Data        : num [1:2162, 1:2162] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:2162] "S1" "S2" "S3" "S4" ...
  .. .. .. .. ..$ : chr [1:2162] "A1" "A2" "A3" "A4" ...
  .. .. ..@ taxa_are_rows: logi FALSE
  .. .. ..$ dim     : int [1:2] 2162 2162
  .. .. ..$ dimnames:List of 2
  .. .. .. ..$ : chr [1:2162] "S1" "S2" "S3" "S4" ...
  .. .. .. ..$ : chr [1:2162] "A1" "A2" "A3" "A4" ...
  ..@ tax_table: NULL
  ..@ sam_data :'data.frame':	2162 obs. of  1 variable:
Formal class 'sample_data' [package "phyloseq"] with 4 slots
  .. .. ..@ .Data    :List of 1
  .. .. .. ..$ : num [1:2162] 1 0 0 0 0 0 0 0 1 0 ...
  .. .. ..@ names    : chr "y"
  .. .. ..@ row.names: chr [1:2162] "S1" "S2" "S3" "S4" ...
  .. .. ..@ .S3Class : chr "data.frame"
  ..@ phy_tree :List of 5
  .. ..$ edge       : int [1:4322, 1:2] 2163 2164 2164 2163 2165 2166 2167 2168 2168 2169 ...
  .. ..$ edge.length: num [1:4322] 0.5545 0.0501 0.0721 0.2573 0.0254 ...
  .. ..$ Nnode      : int 2161
  .. ..$ node.label : chr [1:2161] "" "0.926" "" "'p__Cyanobacteria'" ...
  .. ..$ tip.label  : chr [1:2162] "A1" "A2" "A3" "A4" ...
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
data(obesity)
## maybe str(obesity) ; plot(obesity) ...
}
\keyword{datasets}
