library(sp)     ## for point.in.polygon
library(stringr)## for loadSolid
library(Rglpk)  ## for CanFitTranslation

############################
### Handling of Polygons ###
############################

ConvexHull <- function(points){
  ### gets a n x 2 matrix of 2d points 
  ### returns a m x 2 matrix of those points which form the convex hull, ordered in counterclockwise direction
  
  if (!is.matrix(points)){stop("Expected matrix as input")}
  if (ncol(points)!=2){stop("Input not of required shape")}
  if (nrow(points)<3){stop("Input not of required shape")}
  
  return(points[rev(chull(points[,1],points[,2])),])
}

RotatePoly = function(polygon,alpha){
  ### gets a n x 2 matrix of 2d points
  ### returns a n x 2 matrix of those points rotated by alpha counter-clockwise
  A=RotationMatrix2d(alpha)
  return(polygon %*% t(A))
}

RotateVector = function(vector,alpha){
  ### gets a sequence of length 2 and returns it rotated counter-clockwise by alpha
  A=RotationMatrix2d(alpha)
  return(as.double(A %*% as.matrix(vector,ncol=1)))
}


TranslatePoly <- function(polygon,vector){
  ### translates all points of a polygon by a vector
  
  if (!is.matrix(polygon)){stop("Expected matrix as input")}
  if (ncol(polygon)!=2){stop("Input not of required shape")}
  if (length(vector)!=2){stop("Vector not of length 2")}
  
  
  polygon[,1]=polygon[,1]+vector[1]
  polygon[,2]=polygon[,2]+vector[2]  
  return(polygon)
}

##############################
### Properties of Polygons ###
##############################

Perimeter = function(polygon){
  ### expects the input to be the vertices of a polygon, in the form of a n x 2 matrix
  ### returns the perimeter of a polygon
  
  if (!is.matrix(polygon)){stop("Expected matrix as input")}
  if (ncol(polygon)!=2){stop("Input not of required shape")}
  if (nrow(polygon)<3){stop("Input not of required shape")}
  
  perimeter = 0
  m=nrow(polygon)
  for(i in 1:(m-1)){
    p1 = polygon[i,]
    p2 = polygon[i+1,]
    perimeter = perimeter+sqrt(sum((p1-p2)^2))
  }
  perimeter = perimeter + sqrt(sum((polygon[1,]-polygon[m,])^2))
  return(perimeter)
}

Area <- function(polygon){
  ### expects the input to be the vertices of a polygon, in the form of a n x 2 matrix
  ### returns the area of a polygon 
  
  ### using the "shoelace formula", also known as "Gauss's area formula"
  
  if (!is.matrix(polygon)){stop("Expected matrix as input")}
  if (ncol(polygon)!=2){stop("Input not of required shape")}
  if (nrow(polygon)<3){stop("Input not of required shape")}
  
  n=nrow(polygon)
  
  A=polygon[1,1]*(polygon[2,2]-polygon[n,2])
  
  for (i in 2:(n-1)){
    A=A+polygon[i,1]*(polygon[i+1,2]-polygon[i-1,2])
  }
  
  A=A+polygon[n,1]*(polygon[1,2]-polygon[n-1,2])  
  
  return(abs(A)/2)
}

Diameter <- function(polygon,isPS=FALSE){
  ### expects the input to be the vertices of a polygon, in the form of a n x 2 matrix
  ### returns the length of the longest straight line possible inside a polygon 
  
  if (!is.matrix(polygon)){stop("Expected matrix as input")}
  if (ncol(polygon)!=2){stop("Input not of required shape")}
  if (nrow(polygon)<3){stop("Input not of required shape")}
  
  NrPoints=nrow(polygon) 
  
  if (IsPolygonPointSymmetric(polygon)|isPS){ ## two most distant points will be antipodal
    maxLength=0
    for (i in 1:(NrPoints/2)){
      l=2*sqrt(sum((polygon[i,]^2)))
      maxLength=max(maxLength,l)
    }
  } else {
    maxLength=0
    for (i in 1:(NrPoints-1)){
      for (j in (i+1):NrPoints){
        l=sqrt(sum((polygon[i,]-polygon[j,])^2))
        maxLength=max(maxLength,l)
      }
    }
  }
  return(maxLength)
}


InCircle <- function(polygon){
  ### expects the input to be the vertices of a polygon, in the form of a n x 2 matrix
  ### returns the radius of the largest circle with (0,0) as origin that fits inside the polygon 
  
  if (!is.matrix(polygon)){stop("Expected matrix as input")}
  if (ncol(polygon)!=2){stop("Input not of required shape")}
  if (nrow(polygon)<3){stop("Input not of required shape")}
  
  n=nrow(polygon)
  minRadius=Inf
  for (i in 1:n){
    P1=polygon[i,]
    P2=polygon[i%%n+1,]
    
    P3=ProjectedOrigin(P1,P2)
    minRadius=min(sqrt(sum(P3^2)),minRadius)
  }
  return(minRadius)
}

PolygonInsidePolygon = function(P,Q){
  ### checks, if all points of the polygon P are strictly inside the polygon Q
  ### both P and Q are matrices with 2 columns
  
  if (!is.matrix(P)){stop("Expected matrix as input")}
  if (!is.matrix(Q)){stop("Expected matrix as input")}
  if (ncol(P)!=2){stop("Input not of required shape")}
  if (ncol(Q)!=2){stop("Input not of required shape")}
  if (nrow(P)<3){stop("Input not of required shape")}
  if (nrow(Q)<3){stop("Input not of required shape")}
  
  for(i in 1:nrow(P)){
    if (point.in.polygon(P[i,1],P[i,2],Q[,1],Q[,2])!=1){
      return(FALSE)
    }
  }
  return(TRUE)
}

IsPolygonPointSymmetric <- function(polygon){
  ### expects the input to be the vertices of a polygon, in the form of a n x 2 matrix
  ### checks, if a given polygon is pointsymmetric w.r.t the origin
  
  if (!is.matrix(polygon)){stop("Expected matrix as input")}
  if (ncol(polygon)!=2){stop("Input not of required shape")}
  if (nrow(polygon)<3){stop("Input not of required shape")}
  
  NrPoints=nrow(polygon)
  
  if (NrPoints%%2!=0){return(FALSE)} ### has an odd number of vertices
  
  return(all(polygon[1:(NrPoints/2),]== -polygon[(NrPoints/2+1):NrPoints,]))
}


IsPolyhedronPointSymmetric <- function(polyhedron){
  ### expects the input to be the vertices of a polyhedron, in the form of a n x 3 matrix
  ### checks, if a given polyhedron is pointsymmetric w.r.t the origin
  
  if (!is.matrix(polyhedron)){stop("Expected matrix as input")}
  if (ncol(polyhedron)!=3){stop("Input not of required shape")}
  if (nrow(polyhedron)<3){stop("Input not of required shape")}
  
  NrPoints=nrow(polyhedron)
  
  if (NrPoints%%2!=0){return(FALSE)} ### has an odd number of vertices
  
  for (i in 1:NrPoints){
    AntipodalPointFound=FALSE
    for (j in 1:NrPoints){
      if (all(polyhedron[i,]== -polyhedron[j,])){
        AntipodalPointFound=TRUE
        break
      }
    }
    if (!AntipodalPointFound){return(FALSE)}
  }
  return(TRUE)
}

###########################
### Projection Matrices ###
###########################

RandomSphericalCoordinates <- function(){
  ## returns the pair (theta,phi) of random spherical coordinates, distributed uniformly on the unit sphere
  
  theta=runif(1,0,2*pi)
  phi=acos(runif(1,-1,1))
  return(c(theta,phi))
}

ProjectionMatrix <- function(theta,phi){
  ### returns the projection matrix defined by the parameters
  ### the mapping is orthogonal to X(theta,phi)=(cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi))
  
  A=matrix(nrow=2,ncol=3)
  A[1,]=c(-sin(theta),cos(theta),0)
  A[2,]=c(-cos(theta)*cos(phi),-sin(theta)*cos(phi),sin(phi))
  
  return(A)
}

RandomProjection <- function(){
  ### returns a random 2x3 projection matrix
  theta_phi=RandomSphericalCoordinates()
  return(ProjectionMatrix(theta_phi[1],theta_phi[2]))
}

########################
### Fitting Polygons ###
########################

CanFitTranslation<- function(P,Q){
  ## Uses algorithm from [Chazelle, 1983, Section 3.1] 
  
  ## checks if P can be translated by a vector to fit inside Q,
  ## returns such a vector if possible
  
  if (!is.matrix(P)){stop("Expected matrix as input")}
  if (!is.matrix(Q)){stop("Expected matrix as input")}
  if (ncol(P)!=2){stop("Input not of required shape")}
  if (ncol(Q)!=2){stop("Input not of required shape")}
  if (nrow(P)<3){stop("Input not of required shape")}
  if (nrow(Q)<3){stop("Input not of required shape")}
  
  p=nrow(P)  
  q=nrow(Q)  
  
  A=matrix(nrow=q,ncol=2)
  b=rep(NA,q)
  
  cis=rep(NA,q)
  
  for (i in 1:q){
    w1=Q[i,]
    w2=Q[i%%q+1,]
    w=w2-w1
    
    if (i==1){
      a=RotateVector(w,-pi/2)
      cis[i]=which.max(P%*%a)
    } else {
      ## here maybe improvement possible
      a=RotateVector(w,-pi/2)
      cis[i]=which.max(P%*%a)
    }
  }
  
  for (i in 1:q){
    w1=Q[i,]
    w2=Q[i%%q+1,]
    w=w2-w1
    
    v1=w1+P[1,]-P[cis[i],]
    a=RotateVector(w,pi/2)
    c=a[1]*v1[1]+a[2]*v1[2]
    
    A[i,]=a
    b[i]=c
  }
  
  res=Rglpk_solve_LP(obj=numeric(2),mat=A,dir=rep(">=",q),rhs=b)
  if (res$status==1){
    return(list("fits"=FALSE))
  } else {
    return(list("fits"=TRUE,"vector"=res$solution-P[1,]))
  }
}

CanFitRotation <- function(P,Q,bothPointSymmetric=FALSE){
  ## checks if P can be rotated around the origin to fit inside Q
  ## returns such a angle for P if possible
  
  if (!is.matrix(P)){stop("Expected matrix as input")}
  if (!is.matrix(Q)){stop("Expected matrix as input")}
  if (ncol(P)!=2){stop("Input not of required shape")}
  if (ncol(Q)!=2){stop("Input not of required shape")}
  if (nrow(P)<3){stop("Input not of required shape")}
  if (nrow(Q)<3){stop("Input not of required shape")}
  
  m=nrow(Q)
  n=nrow(P)
  
  Q_EdgesMidPoints=matrix(nrow=m,ncol=2)
  for (i in 1:m){
    Q_EdgesMidPoints[i,]=ProjectedOrigin(Q[i,],Q[i%%m+1,])
  }
  Q_Rs=apply(Q_EdgesMidPoints,1,len)
  Q_Angles=apply(Q_EdgesMidPoints,1,function(r){atan2(r[2],r[1])})
  
  P_Rs=apply(P,1,len)
  P_Angles=apply(P,1,function(r){atan2(r[2],r[1])})
  
  Angles=suppressWarnings(outer(Q_Rs,P_Rs,function(r,s){return(acos(r/s))})) ## some NA's are created
  
  Angles_Plus=outer(Q_Angles,P_Angles,"-")+Angles
  Angles_Minus=outer(Q_Angles,P_Angles,"-")-Angles
  
  criticalAngles=c(Angles_Plus,Angles_Minus)
  criticalAngles=criticalAngles[!is.na(criticalAngles)] ## get rid of all NAs
  
  if (length(criticalAngles)==0){ ## either P always fits in Q or never
    inBetweenAngles=c(0)
  } else {
    
    criticalAngles=criticalAngles%%(2*pi)
    if (bothPointSymmetric){
      criticalAngles=criticalAngles%%pi
      criticalAngles=round(criticalAngles,15)
    }
    
    criticalAngles=unique(criticalAngles)
    criticalAngles=sort(criticalAngles)
    nrCritical=length(criticalAngles)
    
    inBetweenAngles=(head(criticalAngles,-1)+tail(criticalAngles,-1))/2
    
    if (bothPointSymmetric){
      inBetweenAngles=c(inBetweenAngles,(criticalAngles[1]+pi+criticalAngles[nrCritical])/2)
    } else {
      inBetweenAngles=c(inBetweenAngles,mod2pi((criticalAngles[1]+2*pi+criticalAngles[nrCritical])/2))
    }
  }
  
  
  ### demonstration ###
  #plotPoints<- function(Points,col="black",cex=1,type="p",main=""){
  #  maxes = ceiling(max(abs(Points)))
  #  xlim = c(-maxes-1,maxes+1)
  #  ylim = c(-maxes-1,maxes+1)
  #  plot(c(1),c(1),xlim=xlim,ylim=ylim,type="n",xlab="",ylab="",main=main)
  #  points(Points[,1],Points[,2],xlim=xlim,ylim=ylim,type=type,pch=19,col=col,cex=cex,xlab="",ylab="")
  #}
  #
  #for (i in 1:length(inBetweenAngles)){
  #  angle=criticalAngles[i]
  #  plotPoints(RotatePoly(P,angle),cex=0.1,main=paste("critical",angle))
  #  points(rbind(Q,Q),type="l",col="red")
  #  
  #  angle=inBetweenAngles[i]
  #  plotPoints(RotatePoly(P,angle),cex=0.1,main=paste("bewteen",angle))
  #  points(rbind(Q,Q),type="l",col="red")
  #}
  ### \demonstration ###
  
  p=1   ### which point to verify next
  timeSinceLastError=0
  for (angle2Test in inBetweenAngles){
    rotatedP=RotatePoly(P,angle2Test)
    
    while (TRUE){ ### now checks if every point is inside Q
      if (point.in.polygon(rotatedP[p,1],rotatedP[p,2],Q[,1],Q[,2])!=1){
        timeSinceLastError=0
        break
      } else {
        timeSinceLastError=timeSinceLastError+1
        p=p%%n+1
      }
      
      if (timeSinceLastError>n+2){
        ans= list()
        ans[["fits"]] = TRUE
        ans[["angle"]] = angle2Test
        return(ans)
      }
    }
  }
  ans= list()
  ans[["fits"]] = FALSE
  return(ans) 
}

SolveRupert = function(Points, seed = NA, shrinkFaktor = 1,NrProjections=128, wantPlot=TRUE, NrCycles = 1e8){
  ### Checks if a solid has Rupert's property
  ### if possible, it returns a list of the parameters enconding a solution

  wantNinja=FALSE; range=0.0001 ## if true, searches for solutions with very similar projections
  nrDirections=50               ## only used for non-point-symmetric. Similar to area, perimeter, diameter
  
  if (!is.matrix(Points)){stop("Expected Points to be a matrix")}
  if (nrow(Points)<4){stop("Points not of required shape")}
  if (ncol(Points)!=3){stop("Points not of required shape")}
  
  if(!is.na(seed)){
    set.seed(seed)  
  }
  
  isPS = IsPolyhedronPointSymmetric(Points)
  if (isPS){cat("Polyhedron is pointsymmetric\n")}
  if (!isPS){cat("Polyhedron is not pointsymmetric\n")}
  nrPoints=nrow(Points)
  
  solvedData = list()
  
  for (cycle in 1:NrCycles){
    cat("CycleNr:", cycle,"\n")
    
    #### creating variables
    ProjectedPolygons=array(dim=c(NrProjections, nrPoints, 2))
    NrVerticesPolygons=rep(NA,NrProjections)
    
    ThetaPhis=matrix(ncol=2,nrow=NrProjections)
    Alphas   =rep(NA,NrProjections)
    
    if (isPS){nrMetrics=4}               ## area, diameter, incircle, perimeter
    if (!isPS){nrMetrics=3+nrDirections} ## area, diameter, perimeter + Heights
    Metrics=array(dim = c(NrProjections,nrMetrics))
    
    ### computing and storing thetas, phis and alphas (in case of non-PS)
    if (wantNinja){
      baseCoordinates=RandomSphericalCoordinates()
      baseTheta=baseCoordinates[1]
      basePhi=baseCoordinates[2]
      
      ThetaPhis[,1]=baseTheta+runif(NrProjections,-range,range)
      ThetaPhis[,2]=basePhi+runif(NrProjections,-range,range)
      if (!isPS){Alphas=runif(NrProjections,-range,range)}
      
    } else {
      for (i in 1:NrProjections){
        ThetaPhis[i,]=RandomSphericalCoordinates()
        if (!isPS){Alphas[i]=runif(1,0,2*pi)}
      }
    }
    
    ### computing Projections and Hulls and Metrics
    for (i in 1:NrProjections){
      A=ProjectionMatrix(ThetaPhis[i,1],ThetaPhis[i,2])
      ProjectedPoints = Points%*%t(A)
      polygon=ConvexHull(ProjectedPoints)
      if (!isPS){polygon=RotatePoly(polygon,Alphas[i])}
      
      ProjectedPolygons[i,1:nrow(polygon),]=polygon
      NrVerticesPolygons[i]=nrow(polygon)
      
      if (isPS){
        Metrics[i,1]=Area(polygon)
        Metrics[i,2]=Diameter(polygon,isPS=TRUE)
        Metrics[i,3]=InCircle(polygon)
        Metrics[i,4]=Perimeter(polygon)
      } else {
        Metrics[i,1]=Area(polygon)
        Metrics[i,2]=Diameter(polygon)
        Metrics[i,3]=Perimeter(polygon)
        
        for (d in 1:nrDirections){
          rotatedPolygon=RotatePoly(polygon,pi/nrDirections*d)
          Metrics[i,3+d]=max(rotatedPolygon[,1])-min(rotatedPolygon[,1])
        }
      }
    }
    
    for (P_index in 1:NrProjections){
      P=ProjectedPolygons[P_index,1:NrVerticesPolygons[P_index],]
      
      for (Q_index in 1:NrProjections){
        if (Q_index==P_index){next}
        if (any(Metrics[P_index,]>Metrics[Q_index,])){next}
        
        Q=ProjectedPolygons[Q_index,1:NrVerticesPolygons[Q_index],]*shrinkFaktor    
        
        if (isPS){
          result=CanFitRotation(P,Q,TRUE)
          if (result$fits){
            solvedData[["Solution"]] = TRUE
            solvedData[["ThetaPhiP"]] = ThetaPhis[P_index,]
            solvedData[["ThetaPhiQ"]] = ThetaPhis[Q_index,]
            solvedData[["alpha"]] = result$angle
            solvedData[["vector"]] = c(0,0)
            
            if (wantPlot){
              PlotSolution(Points,solvedData)
            }
            
            return(solvedData)
          }
          
        } else {
          result=CanFitTranslation(P,Q)
          if (result$fits){
            
            solvedData[["Solution"]] = TRUE
            solvedData[["ThetaPhiP"]] = ThetaPhis[P_index,]
            solvedData[["ThetaPhiQ"]] = ThetaPhis[Q_index,]
            solvedData[["alpha"]] = Alphas[P_index]-Alphas[Q_index]
            solvedData[["vector"]] = RotateVector(result$vector,-Alphas[Q_index])
            
            if (wantPlot){
              PlotSolution(Points,solvedData)
            }
            return(solvedData)
          }
          
        }
        
      }
    }
  }
  solvedData[["Solution"]] = FALSE
  return(solvedData)
}

PlotSolution <- function(Points,solution){
  ### creates 2 plots of a solution to the Rupert problem
  ### expects "solution" to be the list of paramenters
  
  if (!is.matrix(Points)){stop("Expected Points to be a matrix")}
  if (nrow(Points)<4){stop("Points not of required shape")}
  if (ncol(Points)!=3){stop("Points not of required shape")}
  
  P=Points%*%t(ProjectionMatrix(solution$ThetaPhiP[1],solution$ThetaPhiP[2]))
  P=RotatePoly(P,solution$alpha)
  P=TranslatePoly(P,solution$vector)
  Q=Points%*%t(ProjectionMatrix(solution$ThetaPhiQ[1],solution$ThetaPhiQ[2]))
  
  ####################
  #### first plot ####
  ####################
  
  plot(Q[,1],Q[,2],type="n",xlab="",ylab="",asp=1,bty="n",axes = F)
  for (h in 1:nrow(Edges)){
    Edge=Edges[h,]
    points(P[Edge,1],P[Edge,2],type="l")
    points(Q[Edge,1],Q[Edge,2],type="l",col="red")
  }
  
  #####################
  #### second plot ####
  #####################
  
  ## tries to classify vertices into visable and hidden vertices
  ptsP=ClassifyPoints(P,Edges)
  ptsQ=ClassifyPoints(Q,Edges)
  
  if ((length(ptsP$Hull)==nrow(P))|(length(ptsQ$Hull)==nrow(P))){
    ## all points are on the convex hull -> no special plot is possible with this method
    return()
  }
  
  
  plot(Q[,1],Q[,2],type="n",xlab="",ylab="",asp=1,bty="n",axes = F)
  for (k in 2:1){
    if (k==1){R=P;pts=ptsP;col="gray"}
    if (k==2){R=Q;pts=ptsQ;col="black"}
    
    EdgesVisible=rep(FALSE,nrow(Edges))
    
    ## Edge is Visible, if at least one of its endpoints is in Front:
    for (h in 1:nrow(Edges)){
      Edge=Edges[h,]
      if ((any(Edge[1]==pts$Front))|(any(Edge[2]==pts$Front))){
        EdgesVisible[h]=TRUE
      }
    }
    
    ## Edge is Visible, if both endpoints are on Hull and does not intersect previous visible line:
    for (h in which(!EdgesVisible)){
      Edge=Edges[h,]
      if ((any(Edge[1]==pts$Hull))&(any(Edge[2]==pts$Hull))){
        ## both points on Hull
        doesIntersectAnything =FALSE
        A1=R[Edge[1],]
        A2=R[Edge[2],]
        for (hh in which(EdgesVisible)){
          if (any(Edge==Edges[hh,1])|any(Edge==Edges[hh,2])){next}
          B1=R[Edges[hh,1],]
          B2=R[Edges[hh,2],]
          if (twoLinesIntersect(A1,A2,B1,B2)){doesIntersectAnything=TRUE;break}
        }
        if (!doesIntersectAnything){EdgesVisible[h]=TRUE}
      }
    }
    
    for (h in 1:nrow(Edges)){
      Edge=Edges[h,]
      if (EdgesVisible[h]){
        points(R[Edge,1],R[Edge,2],type="l",col=col,lwd=2)
      } else {
        points(R[Edge,1],R[Edge,2],type="l",col=col,lty=2,lwd=0.5)
      }
    }
  }
}

PlotSolutionBW <- function(Points,solution){
  ### creates 2 plots of a solution to the Rupert problem
  ### expects "solution" to be the list of paramenters
  
  if (!is.matrix(Points)){stop("Expected Points to be a matrix")}
  if (nrow(Points)<4){stop("Points not of required shape")}
  if (ncol(Points)!=3){stop("Points not of required shape")}
  
  P=Points%*%t(ProjectionMatrix(solution$ThetaPhiP[1],solution$ThetaPhiP[2]))
  P=RotatePoly(P,solution$alpha)
  P=TranslatePoly(P,solution$vector)
  Q=Points%*%t(ProjectionMatrix(solution$ThetaPhiQ[1],solution$ThetaPhiQ[2]))
  
  ####################
  #### first plot ####
  ####################
  
  plot(Q[,1],Q[,2],type="n",xlab="",ylab="",asp=1,bty="n",axes = F)
  for (h in 1:nrow(Edges)){
    Edge=Edges[h,]
    points(P[Edge,1],P[Edge,2],type="l")
    points(Q[Edge,1],Q[Edge,2],type="l",col="red")
  }
  
  #####################
  #### second plot ####
  #####################
  
  ## tries to classify vertices into visable and hidden vertices
  ptsP=ClassifyPoints(P,Edges)
  ptsQ=ClassifyPoints(Q,Edges)
  
  if ((length(ptsP$Hull)==nrow(P))|(length(ptsQ$Hull)==nrow(P))){
    ## all points are on the convex hull -> no special plot is possible with this method
    return()
  }
  
  
  plot(Q[,1],Q[,2],type="n",xlab="",ylab="",asp=1,bty="n",axes = F)
  for (k in 2:1){
    if (k==1){R=P;pts=ptsP;col="black";typ=2}
    if (k==2){R=Q;pts=ptsQ;col="black";typ=1}
    
    EdgesVisible=rep(FALSE,nrow(Edges))
    
    ## Edge is Visible, if at least one of its endpoints is in Front:
    for (h in 1:nrow(Edges)){
      Edge=Edges[h,]
      if ((any(Edge[1]==pts$Front))|(any(Edge[2]==pts$Front))){
        EdgesVisible[h]=TRUE
      }
    }
    
    ## Edge is Visible, if both endpoints are on Hull and does not intersect previous visible line:
    for (h in which(!EdgesVisible)){
      Edge=Edges[h,]
      if ((any(Edge[1]==pts$Hull))&(any(Edge[2]==pts$Hull))){
        ## both points on Hull
        doesIntersectAnything =FALSE
        A1=R[Edge[1],]
        A2=R[Edge[2],]
        for (hh in which(EdgesVisible)){
          if (any(Edge==Edges[hh,1])|any(Edge==Edges[hh,2])){next}
          B1=R[Edges[hh,1],]
          B2=R[Edges[hh,2],]
          if (twoLinesIntersect(A1,A2,B1,B2)){doesIntersectAnything=TRUE;break}
        }
        if (!doesIntersectAnything){EdgesVisible[h]=TRUE}
      }
    }
    
    for (h in 1:nrow(Edges)){
      Edge=Edges[h,]
      if (EdgesVisible[h]){
        points(R[Edge,1],R[Edge,2],type="l",col=col,lty=typ,lwd=3)
      } else {
        points(R[Edge,1],R[Edge,2],type="l",col=col,lty=typ,lwd=1)
      }
    }
  }
}

NieuwlandNumber <- function(Points,solution,upperBound=20){
  ### computes the Nieuwland number of a solution to the Rupert Problem 
  
  PHull=RotatePoly(ConvexHull(Points%*%t(ProjectionMatrix(solution$ThetaPhiP[1],solution$ThetaPhiP[2]))),solution$alpha)
  QHull=ConvexHull(Points%*%t(ProjectionMatrix(solution$ThetaPhiQ[1],solution$ThetaPhiQ[2])))
  
  lowerBound=0
  eps=.Machine$double.eps
  
  while(lowerBound<upperBound-eps){
    midPoint=(lowerBound+upperBound)/2
    if (all(point.in.polygon(PHull[,1]*midPoint+solution$vector[1],PHull[,2]*midPoint+solution$vector[2],QHull[,1],QHull[,2])==1)){
      lowerBound=midPoint
    } else {
      upperBound=midPoint
    }
  }
  
  if (upperBound==20){stop("Nieuwland Number Exceeds Upper Bound")}
  return((lowerBound+upperBound)/2)
}

ImproveNieuwland <- function(Points,solution,NrTries=10^6,wantPrint=TRUE, NrTimeSinceLastImprovment = 10000){
  ### this function tries to improve the current Nieuwland number
  
  isPS=IsPolyhedronPointSymmetric(Points)
  
  BestN=NieuwlandNumber(Points,solution)
  
  range=0.02
  timeSinceLastImprovment=0
  NewTry=solution
  
  for (i in 1:NrTries){
    if (i%%10^4==0){if (wantPrint){cat("Progress",round(i/NrTries*100),"%\n")}}
    timeSinceLastImprovment=timeSinceLastImprovment+1
    
    if (timeSinceLastImprovment>=NrTimeSinceLastImprovment){
      range=range/2
      if (wantPrint){cat("now with range",range,"\n")}
      timeSinceLastImprovment=0
    }
    
    NewTry$ThetaPhiP=solution$ThetaPhiP+runif(2,-range,range)
    NewTry$ThetaPhiQ=solution$ThetaPhiQ+runif(2,-range,range)
    NewTry$alpha=solution$alpha+runif(1,-range,range)
    if (isPS){NewTry$vector=c(0,0)}
    if (!isPS){NewTry$vector=solution$vector+runif(2,-range,range)}  
    
    if (NieuwlandNumber(Points,NewTry)>BestN){
      BestN=NieuwlandNumber(Points,NewTry)
      timeSinceLastImprovment=0
      solution=NewTry
      if (wantPrint){cat("New best:",BestN,"\n")}
    }
  }
  
  solution$ThetaPhiP=solution$ThetaPhiP%%(2*pi)
  solution$ThetaPhiQ=solution$ThetaPhiQ%%(2*pi)
  solution$alpha=solution$alpha%%(2*pi)
  
  return(solution)
}

##########################
### Loading Polyhedron ###
##########################

loadSolid = function(n=""){
  ## loads any Platonic, Archimedean, Catalan and Johnson solid
  ## the global variables "Points" and "Edges" are created
  
  Phi = (1+sqrt(5))/2
  if((n == 1)||(tolower(n) == "tetrahedron")){
    Points<<- variations(1,1,1,permutation = "none",signChanges = "even")
    findEdges(sqrt(8))
    cat("Tetrahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 2)||(toupper(n) == "CUBE")){
    Points<<- variations(1,1,1,permutation = "none",signChanges = "all")
    findEdges(2)
    cat("Cube loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 3)||(toupper(n) == "OCTAHEDRON")){
    Points<<- variations(1,0,0,permutation = "all",signChanges = "all")
    Points<<- uniqueRows(Points)
    findEdges(sqrt(2))
    cat("Octahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 4)||(tolower(n) == "dodecahedron")){
    Points1= variations(1,1,1,permutation = "none",signChanges = "all")
    Points2= variations(0,1/Phi,Phi,permutation = "even",signChanges = "all")
    Points<<-rbind(Points1,Points2)
    Points<<-uniqueRows(Points)
    
    findEdges(2*Phi-2)
    cat("Dodecahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 5)||(tolower(n) == "icosahedron")){
    Points<<- variations(0,Phi,1,permutation = "even",signChanges = "all")
    Points<<-uniqueRows(Points)      
    
    findEdges(2)
    cat("Icosahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 6)||(toupper(n) == "TRUNCATED TETRAHEDRON")){
    Points<<- variations(1,1,3,permutation = "all",signChanges = "even")
    Points<<-uniqueRows(Points)
    
    findEdges(sqrt(8))
    cat("Truncated tetrahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 7)||(tolower(n) == "cuboctahedron")){
    Points<<- variations(1,1,0,permutation = "all",signChanges = "all")
    Points<<- uniqueRows(Points)
    
    findEdges(sqrt(2))
    cat("Cuboctahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 8)||(toupper(n) == "TRUNCATED CUBE")){
    Points<<- variations(1,1,sqrt(2)-1,permutation = "all",signChanges = "all")
    Points<<- uniqueRows(Points)
    
    findEdges(2*sqrt(2)-2)
    cat("Truncated Cube loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 9)||(tolower(n) == "truncated octahedron")){
    Points<<- variations(0,1,2,permutation = "all", signChanges = "all")
    Points<<- uniqueRows(Points)
    
    findEdges(sqrt(2))
    cat("Truncated octahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 10)||(tolower(n) == "rhombicuboctahedron")){
    Points<<- variations(1,1,1+sqrt(2),permutation = "all", signChanges = "all")
    Points<<- uniqueRows(Points)
    
    findEdges(2)
    cat("Rhombicuboctahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 11)||(tolower(n) == "truncated cuboctahedron")){
    Points<<- variations(1,1+sqrt(2),1+2*sqrt(2),permutation = "all",signChanges = "all")
    Points<<- uniqueRows(Points)
    
    findEdges(2)
    cat("Truncated cuboctahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 12)||(toupper(n) == "SNUB CUBE")){
    t=1.83928675521416113255185256465328  ##Tribonacci Constant
    Points1= variations(1,1/t,t,permutation = "even",signChanges = "odd")
    Points2= variations(1,1/t,t,permutation = "odd" ,signChanges = "even")
    Points<<-rbind(Points1,Points2)
    
    findEdges(sqrt(2+2/t^2))
    cat("Snub cube loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 13)||(tolower(n) == "icosidodecahedron")){
    Points1 = variations(0,0,Phi,permutation = "all",signChanges = "all")
    Points2 = variations(0.5,Phi/2,Phi^2/2,permutation = "even",signChanges = "all")
    Points<<- rbind(Points1,Points2)
    Points<<- uniqueRows(Points)
    
    findEdges(1)
    cat("Icosidodecahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 14)||(tolower(n) == "truncated dodecahedron")){
    Points1=variations(0,1/Phi,2+Phi,permutation = "even",signChanges = "all")
    Points2=variations(1/Phi,Phi,2*Phi,permutation = "even",signChanges = "all")
    Points3=variations(Phi,2,Phi+1,permutation = "even",signChanges = "all")
    Points<<-rbind(Points1,Points2,Points3)
    Points<<- uniqueRows(Points)
    
    findEdges(sqrt(5)-1)
    cat("Truncated dodecahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 15)||(tolower(n) == "truncated icosahedron")){
    Points1=variations(0,1,3*Phi,permutation = "odd",signChanges = "all")
    Points2=variations(1,2+Phi,2*Phi,permutation = "odd",signChanges = "all")
    Points3=variations(Phi,2,2*Phi+1,permutation = "odd",signChanges = "all")
    Points<<-rbind(Points1,Points2,Points3)
    Points<<- uniqueRows(Points)
    
    findEdges(2)
    cat("Truncated icosahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 16)||(tolower(n) == "rhombicosidodecahedron")){
    Points1=variations(1,1,Phi^3,permutation = "even",signChanges = "all")
    Points2=variations(Phi^2,Phi,2*Phi,permutation = "even",signChanges = "all")
    Points3=variations(2+Phi,0,Phi^2,permutation = "even",signChanges = "all")
    Points<<-rbind(Points1,Points2,Points3)
    Points<<- uniqueRows(Points)
    
    findEdges(2)
    cat("Rhombicosidodecahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 17)||(tolower(n) == "truncated icosidodecahedron")){
    Points1=variations(1/Phi,1/Phi,3+Phi,permutation = "even",signChanges = "all")
    Points2=variations(2/Phi,Phi,1+2*Phi,permutation = "even",signChanges = "all")
    Points3=variations(1/Phi,Phi^2,-1+3*Phi,permutation = "even",signChanges = "all")
    Points4=variations(2*Phi-1,2,2+Phi,permutation = "even",signChanges = "all")
    Points5=variations(Phi,3,2*Phi,permutation = "even",signChanges = "all")
    
    Points<<-rbind(Points1,Points2,Points3,Points4,Points5)
    Points<<- uniqueRows(Points)
    
    findEdges(2*Phi-2)
    cat("Truncated icosidodecahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if((n == 18)||(tolower(n) == "snub dodecahedron")){
    Xi = (Phi/2+1/2*sqrt(Phi-5/27))^(1/3)+(Phi/2-1/2*sqrt(Phi-5/27))^(1/3)
    Beta = Xi*Phi+Phi^2+Phi/Xi
    Alpha = Xi-1/Xi
    
    a = 2*Alpha
    b = 2
    c = 2*Beta
    Points1=variations(a,b,c,permutation = "even",signChanges = "odd")
    
    a = Alpha+Beta/Phi+Phi
    b = -Alpha*Phi+Beta+1/Phi
    c = Alpha/Phi+Beta*Phi-1
    Points2=variations(a,b,c,permutation = "even",signChanges = "odd")
    
    a = Alpha+Beta/Phi-Phi
    b = Alpha*Phi-Beta+1/Phi
    c = Alpha/Phi+Beta*Phi+1 
    Points3=variations(a,b,c,permutation = "even",signChanges = "odd")
    
    a = -Alpha/Phi+Beta*Phi+1
    b = -Alpha+Beta/Phi-Phi
    c = Alpha*Phi+Beta-1/Phi
    Points4=variations(a,b,c,permutation = "even",signChanges = "odd")
    
    a = -Alpha/Phi+Beta*Phi-1
    b = Alpha-Beta/Phi-Phi
    c = Alpha*Phi+Beta+1/Phi
    Points5=variations(a,b,c,permutation = "even",signChanges = "odd")
    
    Points<<-rbind(Points1,Points2,Points3,Points4,Points5)
    findEdges(6.0437380841)
    cat("Snub dodecahedron loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if ((n>=18)&(n<=123)){ ## Catalan or Johnson
    if ((n>=18)&(n<=31)){ ## Catalan
      path=paste0("Catalan\\",list.files("Catalan\\")[n-18],collapse="")
    } else { ##Johnson
      path=paste0("Johnson\\",list.files("Johnson\\")[n-31],collapse="")      
    }
    
    txt=readLines(path)
    
    Cs=c()
    Edges=matrix(nrow=0,ncol=2)
    Vs = matrix(ncol = 3,nrow = 0)
    for(line in txt){
      
      #finding Constants
      if(substr(line,1,1) == "C" & any(substr(line,8,8) == c(".","0","1","2","3","4","5","6","7","8","9")) 
         & any(substr(line,10,10) == c(".","0","1","2","3","4","5","6","7","8","9"))
         & any(substr(line,11,11) == c(".","0","1","2","3","4","5","6","7","8","9"))
         & any(substr(line,12,12) == c(".","0","1","2","3","4","5","6","7","8","9"))){
        
        Cs = c(Cs,as.double(unlist(strsplit(line,"[=]"))[2]))
      } 
      
      #Finding vertices 
      if(substr(line,1,1) == "V" & any(substr(line,2,2) == as.character(0:9))){
        
        l = str_trim(unlist(strsplit(line,"[(,)]")))
        Vs = rbind(Vs,l[2:4])
      }
      
      ##Finding Edges
      if(substr(line,1,1) == "{"){
        l=substr(line,2,nchar(line)-1)
        face=as.numeric(unlist(strsplit(l,",")))
        for (i in 1:length(face)){
          a=c(face[i],face[i%%length(face)+1])+1
          a=sort(a)
          Edges=rbind(Edges,a)
        }
      }
    }
    
    Edges = uniqueRows(Edges)
    Points = matrix(ncol = 3,nrow = nrow(Vs))
    
    for(i in 1:length(Vs)){ ##iterating over entries in matrix Vs
      if (!grepl("C",Vs[i])){
        Points[i]=as.double(Vs[i])
      } else {
        for(j in 1:length(Cs)){
          for(sign in c("-","")){
            if(Vs[i] == paste0(sign,"C",j-1)){
              Points[i] = Cs[j]
              if(sign == "-"){
                Points[i]=Points[i]*(-1)
              }
            }  
          }
        }
      }
    }
    
    if((sum(is.na(Points))==0)==F){
      return(FALSE)
    }
    Edges <<- Edges
    Points <<- Points
    cat(txt[1]," loaded\n Info: V =",nrow(Points),"E =",nrow(Edges),"\n")
  }
  
  if (n==""){
    cat("Platonic Solids:\n")
    cat(" 1 - Tetrahedron \n")
    cat(" 2 - Cube \n")
    cat(" 3 - Octahedron \n")
    cat(" 4 - Dodecahedron \n")
    cat(" 5 - Icosahedron \n\n")
    cat("Archemedian Solids:\n")
    
    cat(" 6 - Truncated tetrahedron       (solved)\n")
    cat(" 7 - Cuboctahedron               (solved)\n")
    cat(" 8 - Truncated cube              (solved)\n")
    cat(" 9 - Truncated octahedron        (solved)\n")
    cat("10 - Rhombicuboctahedron         (solved)\n")
    cat("11 - Truncated cuboctahedron     (solved)\n")
    cat("12 - Snub cube                   (unsolved)\n")
    cat("13 - Icosidodecahedron           (solved)\n")
    cat("14 - Truncated dodecahedron      (solved)\n")
    cat("15 - Truncated icosahedron       (solved)\n")
    cat("16 - Rhombicosidodecahedron      (unsolved)\n")
    cat("17 - Truncated icosidodecahedron (solved)\n")
    cat("18 - Snub dodecahedron           (unsolved)\n\n")
    
    cat("Catalan Solids:\n")
    for (i in 1:length(list.files("Catalan\\"))){
      a=readLines(paste0("Catalan\\",list.files("Catalan\\")[i],collapse=""))
      cat(i+18,"-",a[1],"\n")
    }
    cat("\nJohnson Solids:\n")
    for (i in 1:length(list.files("Johnson\\"))){
      a=readLines(paste0("Johnson\\",list.files("Johnson\\")[i],collapse=""))
      cat(i+31,"-",a[1],"\n")
    }
  }
}

findEdges <- function(sideLengths=NA){
  ## creates the global variable "Edges" containing all the edges of the polyhedron "Points". It therefore uses the side lengths as input
  
  Edges<<- matrix(ncol=2,nrow=0)
  for (sideLength in sideLengths){
    counter=0
    for (i in 1:(nrow(Points)-1)){
      for (j in (i+1):nrow(Points)){
        d=sqrt(sum((Points[i,]-Points[j,])^2))
        if (abs(d-sideLength)<0.0001){
          Edges<<-rbind(Edges,c(i,j))
          counter=counter+1
        }
      }
    }
  }
}

##############
### Extras ###
##############

variations<- function(a,b,c,permutation="none",signChanges="none"){
  ### creates some of the permutations of the triple (a,b,c), with possible sign changes
  ### used in loadSolid()
  
  if (!any(permutation==c("none","even","odd","all"))){stop("permutation = ",permutation," is not allowed")}
  if (!any(signChanges==c("none","even","odd","all"))){stop("signChanges = ",signChanges," is not allowed")}
  
  ## first dealing when only needing to permute:
  if (signChanges=="none"){
    if (permutation=="none"){
      return(matrix(c(a,b,c),nrow=1,ncol=3))
    }
    if (permutation=="even"){
      return(matrix(c(a,b,c,b,c,a,c,a,b),ncol=3,byrow=TRUE))
    }
    if (permutation=="odd"){
      return(matrix(c(a,c,b,c,b,a,b,a,c),ncol=3,byrow=TRUE))
    }
    if (permutation=="all"){
      return(matrix(c(a,b,c,a,c,b,b,a,c,b,c,a,c,a,b,c,b,a),ncol=3,byrow=TRUE))
    }
  }
  
  ## now the other signChanges cases:
  if (signChanges=="all"){
    M=matrix(nrow=0,ncol=3)
    M=rbind(M,variations( a, b, c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations( a, b,-c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations( a,-b, c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations( a,-b,-c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations(-a, b, c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations(-a, b,-c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations(-a,-b, c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations(-a,-b,-c,permutation = permutation,signChanges = "none"))
    return(M)
  }
  
  if (signChanges=="even"){
    M=matrix(nrow=0,ncol=3)
    M=rbind(M,variations( a, b, c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations( a,-b,-c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations(-a, b,-c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations(-a,-b, c,permutation = permutation,signChanges = "none"))
    return(M)
  }
  
  if (signChanges=="odd"){
    M=matrix(nrow=0,ncol=3)
    M=rbind(M,variations(-a, b, c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations( a,-b, c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations( a, b,-c,permutation = permutation,signChanges = "none"))
    M=rbind(M,variations(-a,-b,-c,permutation = permutation,signChanges = "none"))
    return(M)
  }
}

uniqueRows<- function(Matrix){
  ## returns a matrix that only contains the unique rows of the input matrix
  
  if (nrow(Matrix)<1){return(Matrix)}
  
  Rows2Remove=c()
  for (row2Remove in 2:nrow(Matrix)){
    for (previousRows in 1:(row2Remove-1)){
      if (all(Matrix[row2Remove,]==Matrix[previousRows,])){
        Rows2Remove=c(Rows2Remove,row2Remove)
      }
    }
  }
  
  if (length(Rows2Remove)==0){return(Matrix)} ## there is nothing to remove
  
  M=Matrix[-Rows2Remove,]
  return(M)
}

RotationMatrix2d <- function(alpha){
  ### returns the 2x2 rotation matrix, that rotates by an angle alpha counterclockwise
  A=matrix(nrow=2,ncol=2)
  A[,1]=c(cos(alpha),sin(alpha))
  A[,2]=c(-sin(alpha),cos(alpha))
  
  return(A)
}

len<- function(v){
  ## returns the euclidean norm of a vector
  return(sqrt(sum(v^2)))
}

ScalarProduct <- function(v1,v2){
  ## returns the scalar product of two vectors
  return(sum(v1*v2))
}

ProjectedOrigin<- function (P1,P2){
  ### it returns the projection of the point (0,0) onto the line going through two points P1 and P2
  
  t=-ScalarProduct(P1,P2-P1)/ScalarProduct(P2-P1,P2-P1)
  return(P1+t*(P2-P1))
}

twoLinesIntersect <- function(A1,A2,B1,B2){
  ### returns True, if the lines A1A2 and B1B2 intersect
  
  ##Step 1: all Points will be translated, so that A1=(0,0)
  vec=A1
  A1=A1-vec
  A2=A2-vec
  B1=B1-vec
  B2=B2-vec
  
  ## Step 2: all Points will be rotated, so that A2=(0,y) with y>0
  phi=atan2(A2[2],A2[1])
  A1=RotateVector(A1,-phi+pi/2)
  A2=RotateVector(A2,-phi+pi/2)
  B1=RotateVector(B1,-phi+pi/2)
  B2=RotateVector(B2,-phi+pi/2)
  
  ## Step 3: The equation y=kx+d for the line passing B1 and B2 is calculated
  k=(B2[2]-B1[2])/(B2[1]-B1[1])
  d=B1[2]-k*B1[1]
  
  ## Step 4: Check if certain criteria are met
  
  if (B1[1]*B2[1]>0){return(FALSE)} ## B1 and B2 lie on the same side of the line x=0
  
  return((0<d)& (d<A2[2])) ## check if B1B2 crosses the y-axis on the right height
}

ClassifyPoints <- function(P,Edges){
  ## given a projection of a polyhedra, it classifies its vertices to be either
  ## .) on the convex hull
  ## .) in the front -> visible
  ## .) in the back -> not visible
  
  pts<- list()
  pts$Hull=chull(P)
  
  if (length(pts$Hull)==nrow(P)){
    ## all points are on the convex hull -> no special plot is possible with this method
    return(pts)
  }
  pivotP=(1:nrow(P))[-pts$Hull][1]
  pts$Front=c(pivotP)
  AnyChangesInRound=T
  while(AnyChangesInRound){
    AnyChangesInRound=F
    
    for (i in 1:nrow(Edges)){
      for (j in 1:2){
        v1=Edges[i,j]
        v2=Edges[i,3-j]
        ## test if v1 is in pts$Front and v2 can be added
        if (!any(pts$Front==v1)){next} ## v1 not in pts$Front
        if (any(pts$Front==v2)){next}  ## v2 already in pts$Front
        if (any(v2==pts$Hull)){next}     ## v2 on hull
        pts$Front=c(pts$Front,v2)
        AnyChangesInRound=T
      }
    }
  }
  
  pts$Back=(1:nrow(P))[-c(pts$Hull,pts$Front)]
  
  return(pts)
}

mod2pi <- function(v){
  ### computes v mod 2*pi in [0,2pi)
  return(v%%(2*pi))
}