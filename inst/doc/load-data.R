## ----setup---------------------------------------------------------------
baseLoc <- system.file(package="riskSharing")
extPath <- file.path(baseLoc, "extdata")
dataPath <- file.path("data")

## ----read-and-save, eval=FALSE-------------------------------------------
#  nyakatoke <- read.csv(file.path(extPath,"tanzania_data.csv"))
#  save(nyakatoke, file = file.path(dataPath, "nyakatoke.RData"))

