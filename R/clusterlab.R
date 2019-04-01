#' clusterlab
#'
#' This function runs clusterlab which is a simulator for Gaussian clusters. The default method positions cluster
#' centers on the perimeter of a circle, before creating gaussian clusters around them and projecting the 2D
#' co-ordinates into high dimensional feature space. This method allows control over the spacing, variance, and size
#' of the clusters. Also included is a simple random cluster simulator where the spacing of the clusters cannot be 
#' controlled precisely, but the other parameters can.
#' 
#' @param centers Numerical value: the number of clusters to simulate (N)
#' @param r Numerical value: the number of units of the radius of the circle on which the clusters are generated
#' @param sdvec Numerical vector: standard deviation of each cluster, N values are required
#' @param alphas Numerical vector: how many units to push each cluster away from the initial placement, N values are required
#' @param centralcluster Logical flag: whether to place a cluster in the middle of the rest
#' @param numbervec Numerical vector: the number of samples in each cluster, N values are required
#' @param features Numerical value: the number of features for the data
#' @param seed Numerical value: fixes the seed if you want to repeat results, set the seed to 123 for example here
#' @param rings Numerical value: the number of concentric rings to generate (previous settings apply to all ring clusters)
#' @param ringalphas Numerical vector: a vector of numbers to push each ring out by, must equal number of rings
#' @param ringthetas Numerical vector: a vector of angles to rotate each ring by, must equal number of rings
#' @param outliers Numerical value: the number of outliers to create
#' @param outlierdist Numerical value: a distance value to move the outliers by
#' @param mode Character string: whether to use the standard method (circle), or simple random placement (random)
#' @param minalloweddist Numerical value: minimum distance between the randomised cluster centers, otherwise repeat randomisation
#' @param pcafontsize Numerical value: the font size of the pca
#' @param showplots Logical flag: whether to remove the plots
#'
#' @return A list, containing: 
#' 1) the synthetic data
#' 2) cluster membership matrix
#' @export
#'
#' @examples
#' synthetic <- clusterlab(centers=4,r=8,sdvec=c(2.5,2.5,2.5,2.5),   
#' alphas=c(1,1,1,1),centralcluster=FALSE,   
#' numbervec=c(50,50,50,50)) # for a six cluster solution)   
#' 

clusterlab <- function(centers=1,r=8,sdvec=NULL,alphas=NULL,centralcluster=FALSE,
                       numbervec=NULL,features=500,seed=NULL,rings=NULL,ringalphas=NULL,
                       ringthetas=NULL,outliers=NULL,outlierdist=NULL,mode=c('circle','random'),
                       minalloweddist=0,pcafontsize=18,showplots=TRUE){
  
  mode <- match.arg(mode)
  message('***clusterlab***')
  message(paste('mode:',mode))
  message('simulating clusters...')
  
  if (is.numeric(seed)==TRUE){
    set.seed(seed)
  }
  if (is.null(sdvec) == TRUE){
    if (mode == 'circle'){
      message('user has not set standard deviation of clusters, setting automatically...')
      sdvec <- rep(1,centers)
    }else if (mode == 'random'){
      message('user has not set standard deviation of clusters, setting automatically...')
      sdvec <- rep(0.1,centers)
    }
  }
  if (length(numbervec) != centers){
    message('user has not set length of numbervec equal to number of clusters, setting automatically...')
    numbervec <- rep(100,centers)
  }
  if (length(sdvec) != centers){
    if (mode == 'circle'){
      message('user has not set length of sdvec equal to number of clusters, setting automatically...')
      sdvec <- rep(1,centers)
    }else if (mode == 'random'){
      message('user has not set length of sdvec equal to number of clusters, setting automatically...')
      sdvec <- rep(0.1,centers)
    }
  }
  
  if (mode == 'circle'){
    if (centers != 1){
      if (is.null(alphas) == TRUE){ # user does not give extra annotation data
        message('user has not set alphas of clusters, setting automatically...')
        alphas <- rep(1,centers)
      }
    }
    if (abs(max(sdvec) - min(sdvec)) != 0 & is.null(rings) == FALSE){
      message('multiple ring method does not currently allow individual cluster variances, setting equal to first element')
      sdvec <- rep(sdvec[1],centers)
    }
    if (is.null(alphas) == TRUE){
      message('user has not set alphas of clusters, setting automatically...')
      alphas <- rep(1,centers)
    }
    if (length(alphas) != centers){
      message('user has not set length of alphas equal to number of clusters, setting automatically...')
      alphas <- rep(1,centers)
    }
    if (abs(max(alphas) - min(alphas)) != 0 & is.null(rings) == FALSE){
      message('multiple ring method does not currently allow individual alphas, setting equal to first element')
      message('please use the ringalphas parameter to seperate the rings by a specified degree...')
      alphas <- rep(alphas[1],centers)
    }
    if (abs(max(numbervec) - min(numbervec)) != 0 & is.null(rings) == FALSE){
      message('multiple ring method does not currently allow individual cluster sizes, setting equal to first element')
      numbervec <- rep(numbervec[1],centers)
    }
    if (is.null(rings) == FALSE & centralcluster == TRUE){
      message('ring method does not currently allow a central cluster to be generated, skipping')
    }
    if (is.null(rings) == FALSE & is.null(ringalphas) == TRUE){
      message('ring alphas not set, setting automatically...')
      ringalphas <- seq(2,rings*2,2)
    }
    if (is.null(rings) == FALSE & is.null(ringthetas) == TRUE){
      message('ring thetas not set, setting automatically...')
      ringthetas <- rep(0,rings)
    }
    if (is.null(outliers) == FALSE & is.null(outlierdist) == TRUE){
      message('outlier degree not set, setting automatically...')
      outlierdist <- 20
    }
    if (is.null(rings) == TRUE){ # doing without rings 
      if (centers != 1){ # we are generating more than one cluster
        # create N points on a circle in 2D space that are evenly spaced which are samples
        if (centralcluster==TRUE){ # if we are having a cluster at 0,0
          matrix <- matrix(nrow=centers-1,ncol=2)
          n = centers-1
          i = 1
          for (x in seq(0,n-1)){ 
            x1 <- cos(2*pi/n*x)*r
            y1 <- sin(2*pi/n*x)*r
            matrix[i,1] <- x1
            matrix[i,2] <- y1
            i = i + 1
          }
          matrix <- rbind(matrix,c(0,0)) # add a center at the 0,0 co ordinate 
        }else{ # no central cluster
          matrix <- matrix(nrow=centers,ncol=2)
          n = centers
          i = 1
          for (x in seq(0,n-1)){ # edited to N instead of N-1
            x1 <- cos(2*pi/n*x)*r
            y1 <- sin(2*pi/n*x)*r
            matrix[i,1] <- x1
            matrix[i,2] <- y1
            i = i + 1
          }
        }
        # use scalar multiplication to pull the centers further apart by a parameter, alpha
        matrix2 <- matrix(nrow=nrow(matrix),ncol=2)
        for (row in seq(1,nrow(matrix))){
          matrix2[row,] <- matrix[row,]*alphas[row] # scalar multiplication
        }
        # using rnorm create new points based on each center with SD=X,mean=0 (shifting point then saving)
        newsamplematrix <- matrix(ncol=2,nrow=0)
        identitymatrix <- matrix(ncol=2,nrow=sum(numbervec)) # hold the cluster identity (ground truth)
        c <- 1
        for (j in seq(1,nrow(matrix2))){ # for each center
          z <- numbervec[j] # select the number of samples for each cluster
          sd <- sdvec[j]
          for (sample in seq(1,z)){ # loop to create new samples by adding noise
            newsample <- matrix2[j,]+rnorm(2,mean=0,sd=sd) # add the noise to the cluster center
            newsamplematrix <- rbind(newsamplematrix,newsample) # add the newsample to the new sample matrix
            identitymatrix[c,1] <- paste('c',j,'s',sample,sep='')
            identitymatrix[c,2] <- j
            c = c + 1
          }
        }
        identitymatrix <- data.frame(identitymatrix)
        colnames(identitymatrix) <- c('sampleID','cluster')
      }else{ # this is just for one cluster
        matrix2 <- matrix(nrow=1,ncol=2)
        matrix2[1,1] <- 0
        matrix2[1,2] <- 0
        # using rnorm create new points based on each center with SD=X,mean=0 (shifting point then saving)
        newsamplematrix <- matrix(ncol=2,nrow=0)
        identitymatrix <- matrix(ncol=2,nrow=sum(numbervec)) # hold the cluster identity (ground truth)
        c <- 1
        for (j in seq(1,nrow(matrix2))){ # for each center
          z <- numbervec[j] # select the number of samples for each cluster
          sd <- sdvec[j]
          for (sample in seq(1,z)){ # loop to create new samples by adding noise
            newsample <- matrix2[j,]+rnorm(2,mean=0,sd=sd) # add the noise to the cluster center
            newsamplematrix <- rbind(newsamplematrix,newsample) # add the newsample to the new sample matrix
            identitymatrix[c,1] <- paste('c',j,'s',sample,sep='')
            identitymatrix[c,2] <- j
            c = c + 1
          }
        }
        identitymatrix <- data.frame(identitymatrix)
        colnames(identitymatrix) <- c('sampleID','cluster')
      }
    }else{ # otherwise we are generating rings
      message('we are generating clusters arranged in rings...')
      # make a new alpha vector for each ring change alpha
      alphalist <- list()
      alphasoriginal <- alphas
      for (ring in seq(1,rings)){
        #alphas <- c(alphas,alphasoriginal*rings) # choose the ring pull apart factor here
        alphalist[[ring]] <- alphasoriginal*ringalphas[ring]
      }
      ringmatrix <- matrix(nrow=0,ncol=2) # to hold all the co ordinates of all points
      for (ring in seq(1,rings)){ # for every ring create N points
        matrix <- matrix(nrow=centers,ncol=2) # rings does not currently support central cluster
        n = centers
        i = 1
        for (x in seq(0,n-1)){ # edited to N instead of N-1
          x1 <- cos(2*pi/n*x)*r
          y1 <- sin(2*pi/n*x)*r
          matrix[i,1] <- x1
          matrix[i,2] <- y1
          # if we are rotating
          if (is.null(ringthetas)==FALSE){
            matrix[i,1] <- x1*cos(ringthetas[ring])-y1*sin(ringthetas[ring]) # x co-ordinate
            matrix[i,2] <- x1*sin(ringthetas[ring])+y1*cos(ringthetas[ring]) # y co-ordinate
          }
          i = i + 1
        }
        # use scalar multiplication to pull the centers further apart by a parameter, alpha
        # need to create new vector of alphas because of the rings
        alphas <- alphalist[[ring]]
        matrix2 <- matrix(nrow=nrow(matrix),ncol=2)
        for (row in seq(1,nrow(matrix))){
          matrix2[row,] <- matrix[row,]*alphas[row] # scalar multiplication
        }
        ringmatrix <- rbind(ringmatrix, matrix2)
      }
      
      matrix2 <- ringmatrix # over writing matrix2
      numbervec <- rep(numbervec,rings) # increase numbervec because of rings
      sdvec <- rep(sdvec,rings) # increase sdvec because of rings
      
      # using rnorm create new points based on each center with SD=X,mean=0 (shifting point then saving)
      newsamplematrix <- matrix(ncol=2,nrow=0)
      identitymatrix <- matrix(ncol=2,nrow=sum(numbervec)) # hold the cluster identity (ground truth)
      c <- 1
      for (j in seq(1,nrow(matrix2))){ # for each center
        z <- numbervec[j] # select the number of samples for each cluster
        sd <- sdvec[j]
        for (sample in seq(1,z)){ # loop to create new samples by adding noise
          newsample <- matrix2[j,]+rnorm(2,mean=0,sd=sd) # add the noise to the cluster center
          newsamplematrix <- rbind(newsamplematrix,newsample) # add the newsample to the new sample matrix
          identitymatrix[c,1] <- paste('c',j,'s',sample,sep='')
          identitymatrix[c,2] <- j
          c = c + 1
        }
      }
      identitymatrix <- data.frame(identitymatrix)
      colnames(identitymatrix) <- c('sampleID','cluster')
    }
    
    if (is.null(outliers) == FALSE){ 
      message('we are generating outliers...')
      # find the centroid, as we only want to shift the outside points, and bind as the last column
      res <- newsamplematrix # this holds all of the data points for all clusters
      res2 <- t(res)
      test <- cbind(res2,rowMeans(res2)) # find centroid
      test <- data.frame(test)
      colnames(test) <- identitymatrix$sampleID
      # compute distance matrix
      colnames(test)[ncol(test)] <- 'centroid'
      testdf <- reshape::melt(as.matrix(dist(t(test))), varnames = c("row", "col"))
      # get only the distances from the centroid
      testdf <- subset(testdf,testdf$row=='centroid')
      # convert to data frame and label
      mydata <- as.data.frame(res)
      mydata <- data.frame(t(mydata))
      colnames(mydata) <- identitymatrix$sampleID
      # make the outliers
      d <- outlierdist
      for (outlier in seq(1,outliers)){ # starting with the most far away, push them out
        faraway <- as.character(testdf$col[which(testdf$value==sort(testdf$value, decreasing=T)[outlier])]) 
        theta <- sample(1:360, 1) # random angle to rotate them by
        #theta <- 0
        mydata[,faraway][1] = mydata[,faraway][1] + d * cos(theta)
        mydata[,faraway][2]= mydata[,faraway][2] + d * sin(theta)
      }
      newsamplematrix <- t(mydata)
    }else{
      matrix2 <- newsamplematrix
    }
    
    # reverse PCA to generate N dimensional synthetic dataset
    n2 <- features
    x <- rnorm(n2, mean = 0, sd = 0.1) # fixed eigenvector1
    y <- rnorm(n2, mean = 0, sd = 0.1) # fixed eigenvector2
    matrix2 <- newsamplematrix
    res = matrix(nrow = nrow(matrix2), ncol = n2)
    for (i in seq(1,nrow(matrix2),1)){
      a <- matrix2[i,1] # get fixed co ordinate1 for entire row
      b <- matrix2[i,2] # get fixed co ordinate2 for entire row
      xk <- x
      yk <- y
      answer <- a*xk + b*yk
      res[i,] <- answer
    }
    
    mydata <- as.data.frame(res)
    pca1 = prcomp(mydata)
    scores <- data.frame(pca1$x) # PC score matrix
    # plot example
    if (showplots == TRUE){
      p <- ggplot(data = scores, aes_string(x = 'PC1', y = 'PC2', colour = identitymatrix$cluster) ) + geom_point(size=3) +
        theme_bw() + 
        theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.y = element_text(size = pcafontsize, colour = 'black'),
              axis.text.x = element_text(size = pcafontsize, colour = 'black'),
              axis.title.x = element_text(size = pcafontsize),
              axis.title.y = element_text(size = pcafontsize))
      print(p)
    }
    mydata <- data.frame(t(mydata))
    colnames(mydata) <- identitymatrix$sampleID
    
    message('finished.')
    
    newlist <- list('synthetic_data' = mydata, 'identity_matrix' = identitymatrix)
    
    return(newlist)
    
    ###############################################################################################
    
  }else if (mode == 'random'){ # this is for the simple random mode
    
    matrix <- matrix(ncol=10,nrow=0)
    identitymatrix <- matrix(ncol=2,nrow=sum(numbervec))
    # make random centers
    c = 1
    for (sample in seq(1,centers)){
      a <- rnorm(n=features) # scale up factor
      suppressWarnings(matrix <- rbind(matrix,a))
      identitymatrix[c,1] <- paste('c',c,'s1',sep='') # s1 because these are the first samples
      identitymatrix[c,2] <- c
      c = c + 1
    }
    row.names(matrix) <- paste('c',seq(1,centers),'s1',sep='') # clusteri, sample1
    d <- matrix
    stored <- dist(d)
    suppressWarnings(mindist <- min(dist(d)))
    
    xi <- 1
    while (mindist < minalloweddist){
      message('randomising centers again as value/s are lower than the set minimum distance...')
      matrix <- matrix(ncol=features,nrow=0)
      for (sample in seq(1,centers)){
        a <- rnorm(n=features) # scale up factor
        matrix <- rbind(matrix,a)
      }
      row.names(matrix) <- paste('c',seq(1,centers),'s1',sep='') # clusteri, sample1
      d <- matrix
      mindist <- min(dist(d))
      xi = xi + 1
      if (xi == 100){
        stop('minimum distance too high for data')
      }
    }
    stored <- dist(d)
    # make Gaussian clusters
    i = 1
    for (center in seq(1,centers)){ # for each cluster in my sequence of clusters
      centeri <- d[center,] # get the cluster central sample
      for (sample in seq(1,numbervec[center]-1)){ # size of ith cluster minus 1
        newsample <- centeri+rnorm(n=length(centeri),sd=sdvec[center])
        d <- rbind(d, newsample)
        row.names(d)[i+centers] <- paste('c',center,'s',sample+1,sep='')
        identitymatrix[i+centers,1] <- paste('c',center,'s',sample+1,sep='') # s1 because these are the first samples
        identitymatrix[i+centers,2] <- center
        i = i + 1
      }
    } 
    identitymatrix <- data.frame(identitymatrix)
    colnames(identitymatrix) <- c('sampleID','cluster')
    #
    mydata <- as.data.frame(d)
    pca1 = prcomp(mydata)
    scores <- data.frame(pca1$x) # PC score matrix
    # plot example
    if (showplots == TRUE){
      p <- ggplot(data = scores, aes_string(x = 'PC1', y = 'PC2', colour = identitymatrix$cluster) ) + geom_point(size=3) +
        theme_bw() + 
        theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.y = element_text(size = pcafontsize, colour = 'black'),
              axis.text.x = element_text(size = pcafontsize, colour = 'black'),
              axis.title.x = element_text(size = pcafontsize),
              axis.title.y = element_text(size = pcafontsize))
      print(p)
    }
    mydata <- data.frame(t(mydata))
    
    newlist <- list('synthetic_data' = mydata, 'identity_matrix' = identitymatrix)
    
    message('finished.')
    
    return(newlist)
  }
}
