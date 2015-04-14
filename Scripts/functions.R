
#so you need to convert all your x and y based on the grain size.Anyways somehow we end up with a matrix of site by species matrix. 
data = read.csv('./data/cross_comms.csv')
grainSize = 156.25
species.RichnessMatrix = as.matrix(data[which(data$grain==grainSize),-(1:3)])

xyCoords = data[,2:3]
dist.Mat = dist(xyCoords)*sqrt(grainSize)

H <- ceiling(dist.Mat/sqrt(grainSize))*grainSize
hmax <- round((max(dist.Mat/2)/grainSize)*grainSize)
H[H>hmax] <- max(H)
dist.Classes = unique(H)
H <- as.matrix(H)

sites.ByDist = list()
for(i in 1:length(dist.Classes))
    sites.ByDist[[i]] = which(H == dist.Classes[i], arr.ind = TRUE)

cov.Mat = list()
for(i in 1:length(dist.Classes)){
    diff = species.RichnessMatrix[sites.ByDist[[i]][,1],]-species.RichnessMatrix[sites.ByDist[[i]][,2],]
    cov.Mat[[i]] = matrix(0, nrow = ncol(diff), ncol = ncol(diff))
    for(j in 1:nrow(diff)){
        cov.Mat[[i]] = cov.Mat[[i]] + diff[j,]%*%t(diff[j,])
    }
}

C = Reduce('+',cov.Mat)

eigen.ForC<-eigen(C)
eigen.values<-eigen.ForC$values
eigen.vectors<-eigen.ForC$vectors

weighted.Eigenvalue = t(eigen.vectors[1])%*%cov.Mat%*%eigen.vectors[1]
