# tree_vario
# Research project for College of Charleston Quantitative Methods class.

#ok so the code should proceed like so...
#1.) Read in and figure out what grain you want to analyze shit on 
#2.) Apply the grain to the xy coordinates of all your data
#3.) Then apply the data to your function. 

dat = read.csv('./data/cross_comms.csv')
head(dat)
# spatial grains are in units of sq meters
unique(dat$grain)

table(dat$grain)
# compute quadrat total abundances across species
quad_abu = apply(dat[ , -(1:3)], 1, sum)
# compute quadrat average density for each grain
tapply(quad_abu, list(dat$grain), mean)
# in this case the middle grain 156.25 seems like a reasonable choice.


vario <- function(object.mso, object.xy, grain){
    xyLabels <- v(match('x',labels(object.mso)), match('y'), labels(object.mso))
    if('NA' %in% xyLabels)
        xyCoords = object.xy
    else
        xyCoords = object.mso[,c(xyLabels)]

#Then the rest of the calculations will be farily simple. 

X <- matrix(runif(16,0,10), ncol = 4)
cor.Mat <- cor(X)
cor.InvMat<- solve(cor.Mat)

dist.Mat <- matrix(0, ncol = nrow(X), nrow = nrow(X))
for (k in 1:nrow(X)){
    for (i in k:nrow(X)){
        dist.Mat[k,i]= sqrt((X[i,1]-X[k,1])^2+(X[i,2]-X[k,2])^2)
    }
}
#determining the h's is going to be difficult. Like can you determine a range of H's to accept for a list of distance classes?
#So I guess if you use wavelets then you are going to find a set of h's such that they mark the dominating distribution of variance. Then perhaps you say 
dist.Classes = c(0,1,2,3,4,5,6,7,8,9,10)
for (i in 1:length(dist.Classes)-1){
    pairs.byDist = matrix()
    pairs.byDist = which(dist.Mat>dist.Classes[i] & dist.Mat<dist.Classes[i+1], arr.ind = T)
    for (j in 1:nrow(pairs.byDist)){
        covariance.byH[i]= covariance.byH[i] + (X[pairs.byDist[j,1],]-X[pairs.byDist[j,2],])%*%t((X[pairs.byDist[j,1],]-X[pairs.byDist[j,2],])*cor.InvMat)
    }
}

#find out a way to find out how many samples occur within each class and create a vector of these numbers then divide this by the scaler n
#apply a weighted sum to the vector above and each covariance.byH

#ok so pretend totalcovariance.C is our matrix we apply PCA on. 

eigen.ForC<-eigen(C)
eigen.values<-eigen.ForC$values
eigen.vectors<-eigen.ForC$vectors

for (i in 1:length(X)){
    val = eigen.values[i]
    vector = eigen.vectors[i]
    weighted.EigenValue = c()
    for (j in 1:length(covariance.byH))
        weighte.EigenValue[j] = t(vector)%*%covariance.byH[j]%*%vector
}
}