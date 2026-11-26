setwd("C:/Users/Carolina/OneDrive - Universidad de Castilla-La Mancha/Trabajo/Extra/Missings/Package/meantobvs")
devtools::load_all()


install.packages("pdftools")
library(pdftools)

path <- system.file("man", "images/R logo/Rlogo.pdf", package = "meantobvs")

#convert from pdf to png
pdftools::pdf_convert(
  pdf = path,
  format = "png",
  dpi = 500,
  filenames = file.path("./man/images/R logo", "Rlogo.png"),
  pages = 1
)

#build hexagon
install.packages("hexSticker")
library(hexSticker)

path.png <- system.file("man", "images/R logo/Rlogo_cut.png", package = "meantobvs")

sticker(path.png, package = "",
        p_size = 0,
        #image
        s_x = 1,
        s_y = 1.1,
        s_width = 0.67,
        #hexagon
        h_color = "#696969",
        h_fill = "#FFF6EA",

        filename = "man/images/R logo/Rhexagon_raw.png")

sticker(path.png, package = "",
        p_size = 0,
        #image
        s_x = 1,
        s_y = 1.1,
        s_width = 0.67,
        #hexagon
        h_color = "#696969",
        h_fill = "#FFF6EA",

        filename = "man/images/R logo/Rhexagon_raw.pdf")

## personalized qr:

library(pdftools)

path <- system.file("man", "images/R logo/Rlogo_bounded.pdf", package = "meantobvs")

#convert from pdf to png
pdftools::pdf_convert(
  pdf = path,
  format = "png",
  dpi = 500,
  filenames = file.path("./man/images/R logo", "Rlogo.png"),
  pages = 1
)

setwd("./man/R logo")

install.packages("qrcode")
library(qrcode)
install.packages("magick")
library(magick)

url <- "https://github.com/camulro/meantobvs"
qr <- qrcode::qr_code(url, ecl = "H") #error correction = 30%

#save this version of the qr
png("qr_raw.png",
  width = 1000,
  height = 1000,
  bg = "white")

par(mar = c(0, 0, 0, 0))
plot(qr)
dev.off()

#read pngs:
qr_raw <- image_read("qr_raw.png")
# logo <- image_read("Rhexagon.png")
logo <- image_read("Rlogo.png")
#resize:
logo <- image_resize(logo, "210x210")

back <- image_blank(width = 230,
                    height = 230,
                    color = "white")

#centered symbol
logo0 <- image_composite(
  back,
  logo,
  operator = "over",
  gravity = "center"
)

# Superponer el recuadro en el centro del QR
qr_final <- image_composite(
  qr_raw,
  logo0,
  operator = "over",
  gravity = "center"
)

# Guardar resultado
image_write(qr_final, "qrlogoR.png", format = "png")



