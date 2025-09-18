\name{cytokine}
\alias{cytokine}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("cytokine")}
\format{
  The format is:
Formal class 'phyloseq' [package "phyloseq"] with 5 slots
  ..@ otu_table:Formal class 'otu_table' [package "phyloseq"] with 2 slots
  .. .. ..@ .Data        : num [1:190, 1:2323] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:190] "S1" "S2" "S3" "S4" ...
  .. .. .. .. ..$ : chr [1:2323] "A1" "A2" "A3" "A4" ...
  .. .. ..@ taxa_are_rows: logi FALSE
  .. .. ..$ dim     : int [1:2] 190 2323
  .. .. ..$ dimnames:List of 2
  .. .. .. ..$ : chr [1:190] "S1" "S2" "S3" "S4" ...
  .. .. .. ..$ : chr [1:2323] "A1" "A2" "A3" "A4" ...
  ..@ tax_table: NULL
  ..@ sam_data :'data.frame':	190 obs. of  1 variable:
Formal class 'sample_data' [package "phyloseq"] with 4 slots
  .. .. ..@ .Data    :List of 1
  .. .. .. ..$ : num [1:190] 552 526 1130 1130 680 ...
  .. .. ..@ names    : chr "y"
  .. .. ..@ row.names: chr [1:190] "S1" "S2" "S3" "S4" ...
  .. .. ..@ .S3Class : chr "data.frame"
  ..@ phy_tree :List of 5
  .. ..$ edge       : int [1:4510, 1:2] 2324 2325 2326 2326 2325 2327 2328 2329 2330 2331 ...
  .. ..$ edge.length: num [1:4510] 3.49e-02 1.81e-01 5.00e-09 2.37e-03 3.27e-02 ...
  .. ..$ Nnode      : int 2188
  .. ..$ node.label : chr [1:2188] "0.588" "0.985" "0.591" "0.961" ...
  .. ..$ tip.label  : chr [1:2323] "A1" "A2" "A3" "A4" ...
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
data(cytokine)
## maybe str(cytokine) ; plot(cytokine) ...
}
\keyword{datasets}
