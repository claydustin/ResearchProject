source('../Scripts/functions.R')

#setwd('../../data')
data.FileNames = grep("comms", list.files(pattern = "\\.csv$"), value = T)
grainSizes = c(39.06,61.04,39.06,75.56,33.42,39.06,39.06,39.06,87.89,39.06,16.5,61.04,64,56.25,0.25,3.06,39.06,76.56,19.69)
data = list()
for(i in 1:length(data.FileNames))
    data[[i]] = read.csv(data.FileNames[i])

msoList = list()
for(i in 1:length(data))
    msoList[[i]] = mso2(data[[i]][,-c(1:3)], data[[i]][,c(2:3)],grainSizes[i])


#read in data and select data separated at a specific grain size
data = read.csv('./oosting_comms.csv')
grainSize = 256
data = data[which(data$grain==grainSize),]
species.RichnessMatrix = as.matrix(data[,-(1:3)])
abundance.bySpecies = apply(species.RichnessMatrix, 2, sum)#this one
species.RichnessMatrix = species.RichnessMatrix[,which(abundance.bySpecies>5,arr.ind = TRUE)]
species.RichnessMatrix = replace(species.RichnessMatrix, list = which(species.RichnessMatrix > 0,arr.ind = TRUE), 1)#and this one

#global species variance
#Var = sum(cov(species.RichnessMatrix))

#create the distance matrix between sites
xyCoords = data[,2:3]
dist.Mat = dist(xyCoords)*sqrt(grainSize)

#create a list of distance classes 
H <- round(dist.Mat/sqrt(grainSize))*sqrt(grainSize)
hmax <- round((max(dist.Mat/2)/sqrt(grainSize))*sqrt(grainSize))
H[H>hmax] <- max(H)
dist.Classes = unique(H[-which(H == max(H), arr.ind = TRUE)])
H <- as.matrix(H)

#separate site pairs according to distance between them 
sites.ByDist = list()
for(i in 1:length(dist.Classes))
    sites.ByDist[[i]] = which(H == dist.Classes[i], arr.ind = TRUE)

#create a covariance matrix at each distance in H
cov.Mat = list()
for(i in 1:length(dist.Classes)){
    diff = species.RichnessMatrix[sites.ByDist[[i]][,1],]-species.RichnessMatrix[sites.ByDist[[i]][,2],]
    cov.Mat[[i]] = matrix(0, nrow = ncol(diff), ncol = ncol(diff))
    for(j in 1:nrow(diff)){
        cov.Mat[[i]] = cov.Mat[[i]] + diff[j,]%*%t(diff[j,])
    }
    cov.Mat[[i]] = (1/(2*nrow(sites.ByDist[[i]])))*cov.Mat[[i]]
}

#create the overall covariance matrix used in PCA and find the global variance 
C = Reduce('+',cov.Mat)
Var = (1/length(dist.Classes))*sum(C)

#perform PCA analysis
eigen.ForC<-eigen(C)
eigen.values<-eigen.ForC$values
eigen.vectors<-eigen.ForC$vectors

#partition eigenvalue according to amount explained at distance H
weighted.Eigenvalue = c()
for(i in 1:length(cov.Mat)){
    weighted.Eigenvalue[i] = t(eigen.vectors[,1])%*%cov.Mat[[i]]%*%eigen.vectors[,1]
}

#Compute the variance of complementarity at each H. I'm getting the 3rd variance estimate as being comp.Cov -weightedEigenvalues "The observed variance of complementarity without axis 1 appeared tooscillate around its global variance, indicating that all larger scale trend had been accounted for by removing PCA axis 1."
comp.Cov = c()
SR.Cov = c()
PCAremoved.Cov = c()
for(i in 1:length(dist.Classes)){
    comp.Cov[i] = sum(diag(cov.Mat[[i]]))
    SR.Cov[i] = sum(cov.Mat[[i]])
    PCAremoved.Cov[i] = sum(comp.Cov[i])-weighted.Eigenvalue[i]
}
variances = list(comp.Cov, SR.Cov, PCAremoved.Cov)

#Plot the diffence variance measures
xrange = range(dist.Classes)
yrange = range(variances)
plotchar = c(15:17)

plot(xrange, yrange, type = "n", xlab = "Distance", ylab = "Variance")
for (i in 1:3){
    lines(dist.Classes, variances[[i]], type = "b", lwd = 1.5, pch = plotchar[i])
}
abline(h = Var, lty= 4)


legend("bottomright", legend = c("Species richness", "Interspecific associations", "Species composition", "Environmental heterogeneity", "without PCA 1", "Single-species aggregation", "global variance without PCA 1"), pch = c(15, 16, 17), bty = 'n')

