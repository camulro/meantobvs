setwd("C:/Users/Carolina/OneDrive - Universidad de Castilla-La Mancha/Trabajo/Extra/Missings/Package/meantobvs")
devtools::load_all()


# install.packages("pdftools")
library(pdftools)

path <- system.file("images/Rscheme.pdf", package = "meantobvs")

#convert from pdf to png
pdftools::pdf_convert(
  pdf = path,
  format = "png",
  dpi = 500,
  filenames = file.path("./inst/images", "Rscheme.png"),
  pages = 1
)

path <- system.file("images/Rscheme2.pdf", package = "meantobvs")

#convert from pdf to png
pdftools::pdf_convert(
  pdf = path,
  format = "png",
  dpi = 500,
  filenames = file.path("./inst/images", "Rscheme2.png"),
  pages = 1
)

