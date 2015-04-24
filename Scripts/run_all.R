source('./Scripts/functions.R')

setwd('./data')
data.FileNames = grep("comms", list.files(pattern = "\\.csv$"), value = T)
grainSizes = c(39.06,61.04,39.06,75.56,33.42,39.06,39.06,39.06,87.89,39.06,16.5,61.04,64,56.25,0.25,3.06,39.06,76.56,19.69)
data = list()
for(i in 1:2)
    data[[i]] = read.csv(data.FileNames[i])
setwd('../')
msoList = list()
para1 = data[[1]]
para2 = data[[1]][,c(2:3)]
para3 = grainSizes[1]
msoList[[1]] = mso2(para1, para2,para3)

for(i in 1:length(data.FileNames))
   msoList[[i]] = mso2(data[[i]], data[[i]][,c(2:3)], FALSE, grainSizes[i])



