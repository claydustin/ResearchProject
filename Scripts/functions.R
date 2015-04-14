
#so you need to convert all your x and y based on the grain size.Anyways somehow we end up with a matrix of site by species matrix. 
species.RichnessMatrix = 
xyCoords = site.BySpecies$xy
dist.Mat = dist(xyCoords)
sites.ByDist = seperateQuadratsByDistance(dist.Mat, grain)
cov.Mat = list()
for(i in 1:length(sites.ByDist)){
    cov.Mat[i] = sapply(sites.ByDist, computeCovarianceMatrixOfSpeciesRichness(x, species.RichnessMatrix))
}
C = matrixSum(cov.Mat)


#dist.byClass = split(dist.Matrix,H) #OK SO now each observation is put into an associated distance class.

C = matrixSum(cov.MatricesByH)

eigen.ForC<-eigen(C)
eigen.values<-eigen.ForC$values
eigen.vectors<-eigen.ForC$vectors


create a matrix which you continually add to using the innerProd of site pair

##computes the inner prod of two sites separated by distance h. 

matrixSum <- function(x){
    sum = matrix(0, nrow = nrow(x), ncol = ncol(x))
    for(i in 1:length(x)){
        sum = sum + x[i]
    }
    return sum;
}
innerProdOfSitePair <- function(i, j){
    site.Diff = i -j
    cov.matrixForSitePair = t(site.Diff)*site.Diff
    return(cov.matrixForSitePair)
}

seperateQuadratsByDistance(mat, grain){
    
    H <- ceiling(mat/grain)*grain
    hmax <- round((max(mat/2)/grain)*grain)
    H[H>hmax] <- max(H)
    sitesByDist = list()
    dist.Classes = unique(H)
    for(i in 1:length(dist.Classes)){
        sitesByDist[i] = which(H == dist.Classes[i], arr.ind = TRUE)
    }
    return sitesByDist
}

computeCovarianceMatrixOfSpeciesRichness(List, mat){
    sum = matrix()
    for(i in 1:length(List)){
        pair = List[i]
        sum = matrixSum(sapply(List,innerProdOfSitePair(mat[pair[1]], mat[pair[2]])), sum) 
    }
    return sum
}

eigenValueByDistance <- function(cov.Mat, v){
    weighted.Eigenvalue = t(v)%*%cov.Mat%*%v
    return weighte.Eigenvalue
}
