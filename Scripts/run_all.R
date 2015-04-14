#read in data and select data separated at a specific grain size
data = read.csv('./data/cross_comms.csv')
grainSize = 156.25
species.RichnessMatrix = as.matrix(data[which(data$grain==grainSize),-(1:3)])

#create the distance matrix between sites
xyCoords = data[,2:3]
dist.Mat = dist(xyCoords)*sqrt(grainSize)

#create a list of distance classes 
H <- ceiling(dist.Mat/sqrt(grainSize))*grainSize
hmax <- round((max(dist.Mat/2)/grainSize)*grainSize)
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

weighted.Eigenvalue = t(eigen.vectors[1])%*%cov.Mat%*%eigen.vectors[1]

#create matrices of complimentarity (where the covariance matrix is not defined for between species variances)
comp.Cov = list()
for(i in 1:length(cov.Mat)){
    comp.Cov[[i]] = diag(diag(cov.Mat[[i]]))
}


