#read in data and select data separated at a specific grain size
data = read.csv('./data/cross_comms.csv')
grainSize = 156.25
data = data[which(data$grain==grainSize),]
species.RichnessMatrix = as.matrix(data[which(data$grain==grainSize),-(1:3)])

#create the distance matrix between sites
xyCoords = data[,2:3]
dist.Mat = dist(xyCoords)*sqrt(grainSize)

#create a list of distance classes 
H <- ceiling(dist.Mat/sqrt(grainSize))*sqrt(grainSize)
hmax <- round((max(dist.Mat/2)/sqrt(grainSize))*sqrt(grainSize))
H[H>hmax] <- max(H)
dist.Classes = unique(H)
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
}

#create the overall covariance matrix used in PCA
C = Reduce('+',cov.Mat)

#perform PCA analysis
eigen.ForC<-eigen(C)
eigen.values<-eigen.ForC$values
eigen.vectors<-eigen.ForC$vectors

#partition eigenvalue according to amount explained at distance H
weighted.Eigenvalues = list()
for(i in 1:length(cov.Mat)){
    weighted.Eigenvalue = t(eigen.vectors[,1])%*%cov.Mat[[i]]%*%eigen.vectors[,1]
}

#Compute the variance of complementarity at each H
comp.Cov = c()
SR.Cov = c()
PCAremoved.Cov = c()
for(i in 1:length(dist.Classes)){
    comp.Cov[i] = sum(diag(cov.Mat[[i]]))
    SR.Cov[i] = sum(cov.Mat[[i]])
    PCAremoved.Cov[i] = sum(cov.Mat[[i]])-weighted.Eigenvalue[i]
}
variances = list(comp.Cov, SR.Cov, PCAremoved.Cov)

#Plot the difference variance measures
plot(dist.Classes, variances[[1]],pch = 15)
par(new = T)
plot(dist.Classes, variances[[2]], pch = 16)
par(new = T)
plot(dist.Classes, variances[[3]], pch = 17)
legend("bottomright", legend = c("Species richness", "Interspecific associations", "Species composition", "Environmental heterogeneity", "without PCA 1", "Single-species aggregation", "global variance without PCA 1"), pch = c(15, 16, 17), bty = 'n')


