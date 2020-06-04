creation_nvar_SR <- function(W=W){
  # W puede ser una matrix de distancias para establecer el orden
  # W también puede ser un objeto de la clase nb
  
  # profvis({
  repmat = function(X,m,n){
    ## R equivalent of repmat (matlab)
    X<-as.matrix(X)
    mx = dim(X)[1]
    nx = dim(X)[2]
    matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)}
  
  if (class(W)=='knn'){ # & class(W)!='nb' & class(W)!='matrix'
    # W <- knn2nb(W)
    # R <- length(W)
    # NNB <- cbind(1:R,t(matrix(unlist(W),nrow = length(W[[1]]))))
    NNB <- cbind(1:R,W$nn,repmat(-99,dim(W$nn)[1],1))
  }
  if (class(W)=='sf'){
     W <- poly2nb(as(hexs.sf, "Spatial"), queen = FALSE)
     R <- length(W)
     co <- sf::st_coordinates(st_centroid(hexs.sf))
     ang <- useful::cart2pol(co[,1],co[,2],degrees = T)[,2]$theta
     for (kk in 1:R){
     W[[kk]]<-hexs.nb[[kk]][sort(ang[hexs.nb[[kk]]],index.return=TRUE)$ix]
     }
    lnnb <- numeric()
    for (i in 1:R){
      lnnb[i] <- length(W[[i]])
    }
    mlnnb <- max(lnnb)
    NNB <- matrix(0,ncol = (mlnnb+1),nrow = R)
    for (i in 1:R){
      NNB[i,] <- c(i,W[[i]],repmat(-99,1,mlnnb-lnnb[i]))
    }
    NNB <- cbind(NNB,repmat(-99,dim(NNB)[1],1))
  }
  if (class(W)=='matrix'){
      R <- dim(W)[1]
      mlnnb <- max(rowSums(W>0))+1
      NNB <- matrix(0,ncol = (mlnnb+1),nrow = R)
      lnnb <- matrix(0,ncol = 1,nrow = R)
      for (i in 1:R){
        tmp=sort(W[i,],index.return=TRUE)
        a <- (tmp$x>0)*(tmp$ix)
        a <- rev(a[a>0])
        lnnb[i] <- length(a)
        NNB[i,] <- c(i,a,repmat(-99,1,mlnnb-lnnb[i]))
      }  
    }

  B <- NNB
  end <- dim(B)[1]*dim(B)[2]
  B1 <- matrix(t(B),ncol = 1)
  B2 <- cbind(B1[1:(end-1)],B1[2:end])
  B2 <- rbind(B2,c(-99,-99))
  nn <- 0
  
  for (i in 1:dim(B2)[1]){
    # Identifico los que estan en la misma linea
    B3 <- B2
    dB2 <- dim(B)[2]
    k <- floor((i-1)/dB2)+1
    quito1 <- seq(1+(k-1)*dB2,(dB2+(k-1)*dB2))
    # Identifico los que contienen un -99
    quito2 <- ((B3[,2]==-99)*(1:dim(B3)[1]))
    quito2 <- quito2[quito2>0] #quito2(2:end);
    quito3 <- ((B3[,1]==-99)*(1:dim(B3)[1]))
    quito3 <- quito3[quito3>0] #quito2(2:end);
    quito23 <- (c(quito2,quito3))
    quito <- (c(quito1,quito23))
    B3 <- B3[-quito,]
    
    if (sum(quito23==i)==0){
      hh <- rowSums(B3==B2[i,1])+rowSums(B3==B2[i,2])
      nk <- c(sum(hh==0),sum(hh==1),sum(hh==2))
      nn=nn+nk
    }
  }
  return(nn)
  # })
}