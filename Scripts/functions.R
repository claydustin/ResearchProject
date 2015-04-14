
#so you need to convert all your x and y based on the grain size.Anyways somehow we end up with a matrix of site by species matrix. 
a <- ifelse(runif(12,0,1)<.5, 1, 0)
b <- ifelse(runif(12,0,1)<.4, 1, 0)
c <- ifelse(runif(12,0,1)<.6, 1, 0)
d <- ifelse(runif(12,0,1)<.5, 1, 0)
x <- rep.int(c(1,2,3,4),3)
y <- c(1,1,1,1,2,2,2,2,3,3,3,3)
mat= cbind(a,b,c,d,x,y)

species.RichnessMatrix = matrix() 
xyCoords = cbind(mat[,5],mat[,6])
dist.Mat = dist(xyCoords)
grain = 1
sites.ByDist = seperateQuadratsByDistance(dist.Mat, grain)

for(i in 1:length(sites.ByDist)){
    cov.Mat[[i]] = computeCovarianceMatrixOfSpeciesRichness(sites.ByDist[[i]], mat[,-(5:6)])
}

cov.Mat = list()

C = Reduce('+',cov.Mat)

eigen.ForC<-eigen(C)
eigen.values<-eigen.ForC$values
eigen.vectors<-eigen.ForC$vectors

weighted.Eigenvalue = eigenvalueByDistance(C, eigen.vectors[1])


innerProdOfSitePair <- function(i, j){
    site.Diff = i -j
    cov.matrixForSitePair = t(site.Diff)*site.Diff
    return(cov.matrixForSitePair)
}

seperateQuadratsByDistance<-function(mat, grain){
    H <- ceiling(mat/grain)*grain
    hmax <- round((max(mat/2)/grain)*grain)
    #H[H>hmax] <- max(H)
    #sitesByDist = list()
    dist.Classes = unique(H)
    H <- as.matrix(H)
    sitesByDist = list()
    for(i in 1:length(dist.Classes))
        sitesByDist[[i]] = which(H == dist.Classes[i], arr.ind = TRUE)
    return(sitesByDist)
}

computeCovarianceMatrixOfSpeciesRichness<- function(List, mat){
    sum = matrix(0, nrow = ncol(mat), ncol = ncol(mat))
    for(i in 1:nrow(List)){
        pair = List[i,]
        site.Diff = mat[pair[1],]-mat[pair[2],]
        sum = multiplyVectors(site.Diff, site.Diff)+ sum 
    }
    return(sum)
}

eigenvalueByDistance <- function(cov.Mat, v){
    weighted.Eigenvalue = t(v)*cov.Mat*v
    return weighted.Eigenvalue
}

multiplyVectors<-function(x,y){
    mat = matrix(0,nrow = length(x), ncol = length(x))
    for(i in 1:length(x)){
        for (j in 1:length(y)){
            mat[i,j] = x[i]*y[j]
        }
    }
    return(mat)
}
