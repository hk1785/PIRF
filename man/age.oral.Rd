\name{age.oral}
\alias{age.oral}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("age.oral")}
\format{
  The format is:
Formal class 'phyloseq' [package "phyloseq"] with 5 slots
  ..@ otu_table:Formal class 'otu_table' [package "phyloseq"] with 2 slots
  .. .. ..@ .Data        : num [1:217, 1:2319] 0.00 0.00 0.00 2.04e-05 0.00 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:217] "S1" "S2" "S3" "S4" ...
  .. .. .. .. ..$ : chr [1:2319] "A1" "A2" "A3" "A4" ...
  .. .. ..@ taxa_are_rows: logi FALSE
  .. .. ..$ dim     : int [1:2] 217 2319
  .. .. ..$ dimnames:List of 2
  .. .. .. ..$ : chr [1:217] "S1" "S2" "S3" "S4" ...
  .. .. .. ..$ : chr [1:2319] "A1" "A2" "A3" "A4" ...
  ..@ tax_table: NULL
  ..@ sam_data :'data.frame':	217 obs. of  1 variable:
Formal class 'sample_data' [package "phyloseq"] with 4 slots
  .. .. ..@ .Data    :List of 1
  .. .. .. ..$ : int [1:217] 22 28 28 28 26 26 30 30 26 30 ...
  .. .. ..@ names    : chr "y"
  .. .. ..@ row.names: chr [1:217] "S1" "S2" "S3" "S4" ...
  .. .. ..@ .S3Class : chr "data.frame"
  ..@ phy_tree :List of 5
  .. ..$ edge       : int [1:4509, 1:2] 2320 2321 2322 2322 2321 2323 2324 2325 2326 2327 ...
  .. ..$ edge.length: num [1:4509] 3.49e-02 1.81e-01 5.00e-09 2.37e-03 3.27e-02 ...
  .. ..$ Nnode      : int 2191
  .. ..$ node.label : chr [1:2191] "0.588" "0.985" "0.591" "0.961" ...
  .. ..$ tip.label  : chr [1:2319] "A1" "A2" "A3" "A4" ...
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
data(age.oral)
## maybe str(age.oral) ; plot(age.oral) ...
}
\keyword{datasets}
