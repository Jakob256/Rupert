#~~~~~~~~~~~~~~~~~~~~~
## 1. Description ####
#~~~~~~~~~~~~~~~~~~~~~

# This script generates the solution tree for the Noperthedron

# It recursively applies the function "RupertDisprover" to an interval in R^5.
# It determines wether the local or global Theorem can be applied.
# Otherwise, it splits the intervals into different parts and calls itself on
# each of those intervals.

# While this function executes, it stores the nodes and how they can be solved
# in a table "df". However, as accessing elements of a table is very slow,
# we store each column independently. ("T_ID", "T_nodetype",...)

# Finally these columns are combined to the dataframe and can be exported as
# a csv. The code written in Sage can then verify that the solution tree is
# valid.


#~~~~~~~~~~~~~~~~~
## 2. Imports ####
#~~~~~~~~~~~~~~~~~


rm(list = ls())     ## clearing the working space
library(sp)         ## for point.in.polygon
library(data.table) ## for writing the large database
library(gmp)        ## for dealing with large integers

#~~~~~~~~~~~~~~~~~~~~~~~~~
## 3. Small functions ####
#~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 3.1 Dealing with Intervals ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

intLen <- function(interval){
  ## returns the length of an interval:
  return(max(interval)-min(interval))
}

midPoint <- function(interval){
  ## returns the length of an interval
  return(sum(interval)/2) 
}

splitInterval <- function(interval,i,nrParts){
  ## Examples: 
  ## splitInterval(c(2,8),1,3)=c(2,4)
  ## splitInterval(c(2,8),2,3)=c(4,6)
  ## splitInterval(c(2,8),3,3)=c(6,8)
  
  step=(interval[2]-interval[1])/nrParts
  return(interval[1]+(i-1)*step+c(0,step))
}


duration2string <- function(s){
  d=s %/% 86400
  h=(s %% 86400) %/% 3600
  m=(s %% 3600) %/% 60
  s=s %% 60
  sprintf("%d days, %02d:%02d:%02d", d, h, m, s)
}


myRound <- function(n,digits){
  return(substr(as.character(n),1,digits+2))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 3.2 Rational approximations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


is.int <- function(a){
  a=a[!is.na(a)]
  return(all(abs(round(a)-a)<0.0001))
}

approximate_rational <- function(x, precision = 0.00001) {
  ## Given a number x, it returns an integer sequence r of length 2,
  ## such that r[1]/r[2] is approximately x.
  
  ## The absolute error will be less than the parameter "precision".
  
  
  if (x==0){return(c(0,1))}
  initial_x = x
  
  a <- floor(x)
  
  ## initial approximation
  num1= 1
  den1= 0

  ## second approximation  
  num2= floor(x)
  den2= 1
  
  x= x-num2/num1
  
  while(TRUE){
    fraction_approximation = num2/den2
    if (abs(fraction_approximation - initial_x) <= precision){break}
    if (x==0){break}

    x = 1 / x
    a = floor(x)
    num1_old = num1
    den1_old = den1
    num1 = num2
    den1 = den2
    num2 = a * num2 + num1_old
    den2 = a * den2 + den1_old
    x = x - a
  }
  g = gcd(num2, den2)
  return(c(num2/g, den2/g))
}

approx_on_circle_helper = function(P,precision=0.00001){
  ## given a point P in R^2, such that P is close to 
  ## the unit circle, this function tries to find a rational
  ## point that is on the unit circle and very close to P
  
  ## it does this by drawing a line l through (1,0) and P. It
  ## then calculates the intersection of this line with the
  ## x-Axis at a point (t,0). It then finds a rational 
  ## approximation of t called t1/t2.
  
  ## It then calculates the intersection of the line through
  ## (0,1) and the point (0,t1/t2) and the unit circle. One can
  ## verify, that this point has the coordinates
  
  ## (2*t1*t2 , t1^2-t2^2) / (t1^2+t2^2)
  
  x = P[1]
  y = P[2]
  t = x/(1-y)
  
  t_rational = approximate_rational(t, precision=precision)
  t1=t_rational[1]
  t2=t_rational[2]
  
  ans = c(2*t1*t2, t1^2-t2^2, t1^2+t2^2)
  
  if(any(ans>9*10^15)){
    t_rational = as.bigz(approximate_rational(t, precision=precision))
    t1=t_rational[1]
    t2=t_rational[2]
    
    ans = c(2*t1*t2,t1^2-t2^2,t1^2+t2^2)
  }
  g = gcd(gcd(ans[1],ans[2]),ans[3])
  
  ## ... after doing the gcd calculation
  ## we should check if bigz is still necessary! 
  
  return(ans/g) #numer1,numer2,common_denom
}


approx_on_circle <- function(P){
  ## wrapper for the previous function
  precision=1e-9
  while(TRUE){
    res=approx_on_circle_helper(P,precision)
    if (all(abs(res)<10^15)){
      return(as.double(res))
    }
    precision=precision*2
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 3.3 Algebtra and Geometry ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ScalarProduct <- function(v1,v2){return(sum(v1*v2))}

len <- function(v){
  ## the euclidean norm of a vector
  return(sqrt(sum(v^2)))
}

len_matrix <- function(v){
  ## returns the euclidean norm of one or more vectors
  return(sqrt(rowSums(v^2)))
}

Rz <- function(a){
  ## The rotation matrix around the z-axis, as defined in the paper
  return(matrix(c(cos(a),sin(a),0,-sin(a),cos(a),0,0,0,1),3,3))
}

M <- function(theta,phi){
  ### returns the projection matrix defined by the parameters
  ### the mapping is orthogonal to X(theta,phi)=(cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi))
  
  A=matrix(nrow=2,ncol=3)
  A[1,]=c(-sin(theta),cos(theta),0)
  A[2,]=c(-cos(theta)*cos(phi),-sin(theta)*cos(phi),sin(phi))
  
  return(A)
}

X <- function(theta,phi){c(cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi))}

R <- function(alpha){
  ### returns the 2x2 rotation matrix, that rotates by an angle alpha counterclockwise
  A=matrix(nrow=2,ncol=2)
  A[1,]=c(cos(alpha),-sin(alpha))
  A[2,]=c(sin(alpha),cos(alpha))
  
  return(A)
}
M_theta_prime <- function(theta,phi,alpha){
  A=matrix(nrow=2,ncol=3)
  
  A[1,]=c(-cos(theta),-sin(theta),0)
  A[2,]=c(sin(theta)*cos(phi),-cos(theta)*cos(phi),0)
  
  return(A)
}


M_phi_prime <- function(theta,phi){
  A=matrix(nrow=2,ncol=3)
  
  A[1,]=c(0,0,0)
  A[2,]=c(cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi))
  
  return(A)
}

R_alpha_prime <- function(alpha){
  A=matrix(nrow=2,ncol=2)
  
  A[,1]=c(-sin(alpha),cos(alpha))
  A[,2]=c(-cos(alpha),-sin(alpha))
  return(A)  
}


ConvexHull <- function(points){
  ### gets a n x 2 matrix of 2d points 
  ### returns a m x 2 matrix of those points which form the convex hull, ordered in counterclockwise direction
  return(points[rev(chull(points[,1],points[,2])),])
}

#~~~~~~~~~~~~~~~~~
## 4. Oracles ####
#~~~~~~~~~~~~~~~~~


masterAnalyserWrapper <- function(Points,T1,V1,T2,V2,A){
  t1=midPoint(T1)
  v1=midPoint(V1)
  t2=midPoint(T2)
  v2=midPoint(V2)
  a=midPoint(A)
  
  P=Points%*%t(M(t1,v1))
  P=P%*%t(R(a))
  Q=Points%*%t(M(t2,v2))
  
  QHull=ConvexHull(Q)
  PHull=ConvexHull(P)
  
  Outside=PHull[point.in.polygon(PHull[,1],PHull[,2],
                                       QHull[,1],QHull[,2])!=1,]
  
  Outside=matrix(Outside,ncol=2) ##I am scared:(
  if (nrow(Outside)==0){master<<-matrix(nrow=0,ncol=2);return()}

  master<<-masterAnalyser(Outside,QHull)
  master<<-matrix(master,ncol=2) ##I am scared:(
  master<<-master[master[,1]>=0,] ##only keep half
  master<<-matrix(master,ncol=2)
  ## warning!!! I used "as.matrix" instead of "matrix" and it broke everything, see "as.matrix(c(1,2),nrow=1)"
  if (nrow(master)==1){return()}
  master<<-master[order(len_matrix(master),decreasing = T),]
}


masterAnalyser <- function(Outside,P){
  ## This function is given a set of points "Outside" and a second set "P".
  ## All the points of "Outside" are outside of "P".
  
  ## For each of the outside points, we calcalate the minimal distance to P.
  ## Specifically, the minimal vectors to P will be returned.
  
  ## Hence this function returns a nrow(Outside) x 2 matrix
  
  ## This function is extremely efficiently
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Distances to Vertices ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  res1=matrix(nrow=nrow(Outside),ncol=2)
  
  ## each row should be a Point Outside
  ## each column a vertex of P
  vec_x=outer(Outside[,1],P[,1],"-")
  vec_y=outer(Outside[,2],P[,2],"-")
  dists2=vec_x^2+vec_y^2
  
  whichVertexClosest=apply(dists2,1,which.min)
  
  res1[,1]=vec_x[cbind(1:nrow(Outside),whichVertexClosest)] ## wow, this is great R syntax
  res1[,2]=vec_y[cbind(1:nrow(Outside),whichVertexClosest)]
  res1=-res1
  
  
  if (F){ ## a great demonstration
    l=max(abs(Outside),abs(P))
    plot(rbind(P,P),type="l",xlim=c(-l,l),ylim=c(-l,l))
    for (i in 1:nrow(Outside)){
      a=Outside[i,]
      b=res1[i,]
      points(c(a[1],a[1]+b[1]),c(a[2],a[2]+b[2]),type="l",col="red")
    }
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Distances to Edges ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~
  
  res2=matrix(nrow=nrow(Outside),ncol=2)
  
  Edges=matrix(nrow=nrow(P),ncol=4)
  Edges[,1:2]=P
  Edges[,3:4]=P[c(2:nrow(P),1),]
  
  NormalVectors=matrix(nrow=nrow(Edges),ncol=2)
  NormalVectors[,1]=Edges[,2]-Edges[,4]
  NormalVectors[,2]=Edges[,3]-Edges[,1]
  NormalVectors=NormalVectors/len_matrix(NormalVectors)
  
  ## scalar product normals to outsiders
  scalar=outer(Outside[,1],NormalVectors[,1],"*")+outer(Outside[,2],NormalVectors[,2],"*")
  
  ## scalar product normals to edges  
  scalar2=Edges[,1]*NormalVectors[,1]+Edges[,2]*NormalVectors[,2]
  
  ## compute difference of scalar products
  ## this is actually the length of the projections
  scalarDif=t(t(scalar)-scalar2)
  
  ## compute the projections of the Outsiders onto the Edges: 
  proj_x=matrix(rep(Outside[,1],nrow(Edges)),nrow=nrow(Outside))-
    scalarDif*matrix(rep(NormalVectors[,1],each=nrow(Outside)),ncol=nrow(Edges))
  proj_y=matrix(rep(Outside[,2],nrow(Edges)),nrow=nrow(Outside))-
    scalarDif*matrix(rep(NormalVectors[,2],each=nrow(Outside)),ncol=nrow(Edges))
  
  
  ## now, we need to decide, if a projection of the outside points lands on the edges
  ## only then the dists are valid
  ## for that, we calculate if the projection is in between the ends of the edges
  
  EdgesStartX=matrix(rep(Edges[,1],nrow(Outside)),ncol=nrow(Edges),byrow = T)
  EdgesStartY=matrix(rep(Edges[,2],nrow(Outside)),ncol=nrow(Edges),byrow = T)
  EdgesEndX=matrix(rep(Edges[,3],nrow(Outside)),ncol=nrow(Edges),byrow = T)
  EdgesEndY=matrix(rep(Edges[,4],nrow(Outside)),ncol=nrow(Edges),byrow = T)
  
  ## I am not so happy about the following lines
  ## There were cases, where ">" did not suffice and ">=" was needed.
  ## This is the case, when the edge is (almost) perfectly vertical or horizontal
  ## maybe I should include a safety margin...
  ## I am not happy:(
  proj_x_ok=(EdgesStartX<=proj_x & proj_x<=EdgesEndX) | (EdgesStartX>=proj_x & proj_x>=EdgesEndX)
  proj_y_ok=(EdgesStartY<=proj_y & proj_y<=EdgesEndY) | (EdgesStartY>=proj_y & proj_y>=EdgesEndY)
  
  proj_ok=proj_x_ok|proj_y_ok ## after many thoughts, I changed it from "&" to "|"
  ##This should fix the vertical/horizontal problem
  
  ## calculate Valid Distances
  Lengths=abs(scalarDif)
  Lengths[!proj_ok]=Inf
  
  bestEdge=apply(Lengths,1,which.min)
  bestEdgeValid=(apply(Lengths,1,min)!=Inf)
  
  res2[,1]=proj_x[cbind(1:nrow(proj_x),bestEdge)]-Outside[,1]
  res2[,2]=proj_y[cbind(1:nrow(proj_y),bestEdge)]-Outside[,2]
  res2[!bestEdgeValid,]=Inf
  
  
  if (F){ ## a great demonstration
    for (i in 1:nrow(Outside)){
      a=Outside[i,]
      b=res2[i,]
      points(c(a[1],a[1]+b[1]),c(a[2],a[2]+b[2]),type="l",col="red")
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## outputting the results ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  takeEdges=(len_matrix(res1)>len_matrix(res2))
  
  res=res1
  res[takeEdges,]=res2[takeEdges,]
  return(res)
}




strong_Oracle <- function(Points,T1,V1,T2,V2,A){
  # returns indices for P1,P2,P3,Q1,Q2,Q3; the signs s_p and s_q; and round(1000*r)
  # if not found, return 0
  # note: round(1000*r) is returned so that all return values are integers
  
  t1=midPoint(T1)
  v1=midPoint(V1)
  t2=midPoint(T2)
  v2=midPoint(V2)
  a=midPoint(A)
  eps=max(intLen(T1),intLen(V1),intLen(T2),intLen(V2),intLen(A))/2
  
  
  rho=max(len_matrix(Points))
  
  g=2*sqrt(2)*rho^2*eps +2*rho^2*eps^2
  
  M1=M(t1,v1)
  M2=M(t2,v2)
  P=Points%*%t(M1)%*%t(R(a))
  Q=Points%*%t(M2)
  
  X1=X(t1,v1)
  X2=X(t2,v2)
  
  
  ## computing candidates
  PHull_cand=sort(chull(P[,1],P[,2]))
  QHull_cand=sort(chull(Q[,1],Q[,2]))
  
  r=(min(len_matrix(Q[QHull_cand,]))-1.42*rho*eps)/rho*0.99
  r=floor(r*1000)/1000
  
  ## without loss of generality, I want the Points of P to be "in front". i.e.s_P=1
  margin=10^(-7)
  PHull_cand_hv=Points%*%X1>1.42*rho*eps+margin ##
  QHull_cand_hv=abs(Points%*%X2)>1.42*rho*eps+margin ## for now, don't discriminate
  
  
  PHull_cand=intersect(PHull_cand,which(PHull_cand_hv))
  QHull_cand=intersect(QHull_cand,which(QHull_cand_hv))
  
  #QHull_cand_hv=sqrt(Q[,1]^2+Q[,2]^2)>=rho*r+1.42*rho*eps #this check is no longer necessary
  #QHull_cand=intersect(QHull_cand,which(QHull_cand_hv))
  
  if (length(PHull_cand)<3){return(0)}
  if (length(QHull_cand)<3){return(0)}
  
  ## column 1 will have the P-indices
  Pairs=matrix(c(PHull_cand,rep(0,length(PHull_cand))),ncol=2)
  #colnames(Pairs)=c("P_index","Q_index")
  
  
  ## can be replaced by these lines
  dif_x=outer(P[Pairs[,1],1],Q[QHull_cand,1],"-")
  dif_y=outer(P[Pairs[,1],2],Q[QHull_cand,2],"-")
  dists2=dif_x^2+dif_y^2
  Pairs[,2]=QHull_cand[apply(dists2,1,which.min)]
  
  
  
  ## Note from 22.09.2024: The following line does nothing, correct?
  Pairs=Pairs[Pairs[,2]%in%QHull_cand,]
  if (length(Pairs)<3*3){return(0)} ## this is because R is stupid
  
  Pairs=cbind(Pairs,delta=0,type=0,maxDelta=-Inf)
  
  ## delta based on the initial definition:
  Pairs[,3]=len_matrix(P[Pairs[,1],]-Q[Pairs[,2],])/2
  
  ## in Q is in the front or the back:
  Pairs[,4]=-1+2*(Points[Pairs[,2],]%*%X2>0)

  for (i in 1:nrow(Pairs)){
    q_ind=Pairs[i,2]
    
    maxDelta=Inf
    ## pruning:
    if (maxDelta<Pairs[i,3]){break} ## this is from 26.04.2024
    
    M2_asd=M(t2,v2)
    Q=Points%*%t(M2_asd)
    
    tries=1:nrow(Points)
    tries=tries[tries!=q_ind]
    noms=len(Q[q_ind,])^2-Q[tries,]%*%Q[q_ind,]-len_matrix(t(t(Points[tries,])-Points[q_ind,]))  *  (2*sqrt(2)*rho*eps+2*rho*eps^2)
    noms=noms-10000*(10)^(-10)## kappa adjustment
    
    dens=(len(Q[q_ind,])+1.42*rho*eps)*
      (len_matrix(Q[rep(q_ind,length(tries)),]-Q[tries,])+2.84*rho*eps)
    fracs=noms/dens*0.95 ## new kappa adjustment
    
    ## maxDelta=min(maxDelta,(fracs-4.5*eps)*(2*r*rho)/2) ## this is the old version, I suspect it was wrong (today 24.4.2024)
    maxDelta=min(maxDelta,(fracs-4.5*eps/(2*r))*(r*rho))
    
    
    
    #print(fracs)
    
    Pairs[i,5]=maxDelta*0.95 ## kappa adjustment
    
  }
  
  Pairs=Pairs[Pairs[,5]>=Pairs[,3],]
  
  #print(Pairs)
  ## Checkpoint: 10.4
  
  if (length(Pairs)<3*3){return(0)} ## this is because R is stupid
  
  Q=Points%*%t(M2)
  mnk=nrow(Pairs)
  
  for (i in 1:(mnk-1)){
    p1=Pairs[i,1]
    for (j in (i+1):mnk){
      p2=Pairs[j,1]
      
      if (ScalarProduct(R(pi/2)%*%P[p1,],P[p2,])<= g+margin){next}
      
      for (k in (i+1):mnk){
        if (length(unique(Pairs[c(i,j,k),4]))!=1){next}
        
        p3=Pairs[k,1]
        
        if (ScalarProduct(R(pi/2)%*%P[p2,],P[p3,])<= g+margin){next}
        if (ScalarProduct(R(pi/2)%*%P[p3,],P[p1,])<= g+margin){next}
        
        ## make the Q search:
        q1=Pairs[i,2]
        q2=Pairs[j,2]
        q3=Pairs[k,2]
        
        
        if (ScalarProduct(R(pi/2)%*%Q[q1,],Q[q2,])<= g+margin){next}
        if (ScalarProduct(R(pi/2)%*%Q[q2,],Q[q3,])<= g+margin){next}
        if (ScalarProduct(R(pi/2)%*%Q[q3,],Q[q1,])<= g+margin){next}
        
        
        mm=Pairs[c(i,j,k),]
        
        
        if (max(mm[,3])>min(mm[,5])){next}
        
        ## I think we are all clear:
        
        ## Note from 22.09.2024:
        ## turns out, we are not all clear:
        ## this is because sometimes the Pairs-matching above yields strange results
        ## I could be more exclusive when finding those pairs
        ## However, for now I simply check here that L condition
        M1=Points[mm[,1],]
        M2=Points[mm[,2],]
        L=solve(M1,M2)
        if (any(abs(L%*%t(L)-diag(3))>0.00001)){next}
        
        return(c(mm[,1],mm[,2],1,mm[1,4],round(1000*r)))
        
        
      }
    }
  }
  
  return(0)
  
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 5. Local and Global Theorem ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


weakTheorem <- function(Points,T1,V1,T2,V2,A,NNN=1){
  
  if (nrow(master)==0){return(FALSE)}
  vec=master[min(NNN,nrow(master)),]
  
  if(all(vec==0)){return(FALSE)}
  vec=vec/len(vec)
  if (round(len(vec),6)!=1){return(FALSE)}
  
  eps=max(intLen(T1),intLen(V1),intLen(T2),intLen(V2),intLen(A))/2
  
  a=midPoint(A)
  t1=midPoint(T1)
  v1=midPoint(V1)
  t2=midPoint(T2)
  v2=midPoint(V2)
  
  hv=Points%*%t(R(a)%*%M(t1,v1))
  S=Points[which.max(hv%*%vec),]
  
  global <<- c(vec,which.max(hv%*%vec))
  
  #~~~~~~~~~~~~~~~
  ## lets gooo ###
  #~~~~~~~~~~~~~~~
  
  G=ScalarProduct(R(a)%*%M(t1,v1)%*%S,vec)-
    eps*abs(ScalarProduct(R_alpha_prime(a)%*%M(t1,v1)%*%S,vec))-
    eps*abs(ScalarProduct(R(a)%*%M_theta_prime(t1,v1)%*%S,vec))-
    eps*abs(ScalarProduct(R(a)%*%M_phi_prime(t1,v1)%*%S,vec))-
    4.5*eps^2
  
  
  hv1=Points%*%t(M(t2,v2))%*%vec
  hv2=Points%*%t(M_theta_prime(t2,v2))%*%vec
  hv3=Points%*%t(M_phi_prime(t2,v2))%*%vec
  
  hv=hv1+eps*abs(hv2)+eps*abs(hv3)+2*eps^2
  H=max(hv)
  
  margin=10^(-7) ## brutally using one margin
  if (G>H+margin){return(TRUE)}
  return(FALSE)
}

strongTheorem <- function(Points,T1,V1,T2,V2,A){
  # consults strongOracle
  # returns false if not applicable
  # returns true if applicable
  
  ## this is the latest version without the |Qi-A|<2R estimate
  
  res=strong_Oracle(Points,T1,V1,T2,V2,A)
  
  if (res[1]==0){return(FALSE)}
  
  global<<-res
  
  P1_index=res[1]
  P2_index=res[2]
  P3_index=res[3]
  Q1_index=res[4]
  Q2_index=res[5]
  Q3_index=res[6]
  s_p=res[7]
  s_q=res[8]
  r=res[9]/1000
  
  stopifnot(s_p%in%c(-1,1))
  stopifnot(s_q%in%c(-1,1))
  stopifnot((0<r) & (r<1))
  
  P1=Points[P1_index,]
  P2=Points[P2_index,]
  P3=Points[P3_index,]
  Q1=Points[Q1_index,]
  Q2=Points[Q2_index,]
  Q3=Points[Q3_index,]
  
  
  t1=midPoint(T1)
  v1=midPoint(V1)
  t2=midPoint(T2)
  v2=midPoint(V2)
  a=midPoint(A)
  
  eps=max(intLen(T1),intLen(V1),intLen(T2),intLen(V2),intLen(A))/2
  
  X_1=X(t1,v1)
  X_2=X(t2,v2)
  
  rho=max(len_matrix(Points))
  
  
  r1=(R(a)%*%M(t1,v1)%*%P1 -M(t2,v2)%*%Q1)/2
  r2=(R(a)%*%M(t1,v1)%*%P2 -M(t2,v2)%*%Q2)/2
  r3=(R(a)%*%M(t1,v1)%*%P3 -M(t2,v2)%*%Q3)/2
  
  delta=max(len(r1),len(r2),len(r3))
  
  g=2*sqrt(2)*rho^2*eps+2*rho^2*eps^2
  
  #~~~~~~~~~~~~~~~~~~~~~~~  
  ## That L condition ####
  #~~~~~~~~~~~~~~~~~~~~~~~  
  
  M1=rbind(P1,P2,P3)
  M2=rbind(Q1,Q2,Q3)
  
  L=solve(M1,M2)
  
  if (any(abs(L%*%t(L)-diag(3))>0.00001)){return(FALSE)}
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~
  ## That s_p condition ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (s_p*ScalarProduct(X_1,P1)<=sqrt(2)*rho*eps){return(FALSE)}
  if (s_p*ScalarProduct(X_1,P2)<=sqrt(2)*rho*eps){return(FALSE)}
  if (s_p*ScalarProduct(X_1,P3)<=sqrt(2)*rho*eps){return(FALSE)}
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~  
  ## That s_q condition ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~  
  
  if (s_q*ScalarProduct(X_2,Q1)<=sqrt(2)*rho*eps){return(FALSE)}
  if (s_q*ScalarProduct(X_2,Q2)<=sqrt(2)*rho*eps){return(FALSE)}
  if (s_q*ScalarProduct(X_2,Q3)<=sqrt(2)*rho*eps){return(FALSE)}
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  ## Those many inequalities ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (ScalarProduct(R(pi/2)%*%M(t1,v1)%*%P1,M(t1,v1)%*%P2)<=g){return(FALSE)}
  if (ScalarProduct(R(pi/2)%*%M(t1,v1)%*%P2,M(t1,v1)%*%P3)<=g){return(FALSE)}
  if (ScalarProduct(R(pi/2)%*%M(t1,v1)%*%P3,M(t1,v1)%*%P1)<=g){return(FALSE)}
  
  if (ScalarProduct(R(pi/2)%*%M(t2,v2)%*%Q1,M(t2,v2)%*%Q2)<=g){return(FALSE)}
  if (ScalarProduct(R(pi/2)%*%M(t2,v2)%*%Q2,M(t2,v2)%*%Q3)<=g){return(FALSE)}
  if (ScalarProduct(R(pi/2)%*%M(t2,v2)%*%Q3,M(t2,v2)%*%Q1)<=g){return(FALSE)}
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Points are far from the origin ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (len(M(t2,v2)%*%Q1)<=r*rho+1.42*rho*eps){return(FALSE)}
  if (len(M(t2,v2)%*%Q2)<=r*rho+1.42*rho*eps){return(FALSE)}
  if (len(M(t2,v2)%*%Q3)<=r*rho+1.42*rho*eps){return(FALSE)}
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## That rational inequality ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  for (j in 1:3){
    if (j==1){Qi=Q1;As=Points[-Q1_index,]}
    if (j==2){Qi=Q2;As=Points[-Q2_index,]}
    if (j==3){Qi=Q3;As=Points[-Q3_index,]}
    
    ## "As" is now the collection of the 89 remaining points:
    
    noms=len(M(t2,v2)%*%Qi)^2- (As%*%t(M(t2,v2)))%*%(M(t2,v2)%*%Qi)-  
      rho*len_matrix(t(Qi-t(As)))*(2*sqrt(2)*eps+2*eps^2)
    
    dens=(len(M(t2,v2)%*%Qi)+1.42*rho*eps) * (len_matrix(t(t(As)-Qi)%*%t(M(t2,v2)))+2.84*rho*eps)
    fracs=noms/dens
    if (any(fracs<(4.5*rho*eps+2*delta)/(2*r*rho))){stop("THIS SHOULD NOT HAPPEN!!!");return(FALSE)}
  
  }
  
  
  return(TRUE)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~
## 6. Rupert Disprover ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~



RupertDisprover <- function(Points,T1,V1,T2,V2,A,depth=0,ID=1){
  ## returns true if a solution can be disproven
  ## returns false otherwise (or throws an error)
  if (depth==0){
    START<<-as.double(Sys.time())
    CALLED<<-0
    MAXDEPTH<<-0
    done<<-0
    todo<<-intLen(T1)*intLen(V1)*intLen(T2)*intLen(V2)*intLen(A)
    WEAK<<-0
    STRONG<<-0
  }
  
  recovery <<- c(T1,V1,T2,V2,A)
  CALLED<<-CALLED+1
  
  MAXDEPTH<<-max(MAXDEPTH,depth)
  percentage=intLen(T1)*intLen(V1)*intLen(T2)*intLen(V2)*intLen(A)/todo
  
  
  if (CALLED%%100==0){
    passed=round(as.double(Sys.time())-START)
    total=round(passed/done/3600)
    
    cat("\r")
    str=rep("",4)
    str[1]=paste("T1: ",myRound(midPoint(T1),4),";  ",
        "V1: ",myRound(midPoint(V1),4),";  ",
        "T2: ",myRound(midPoint(T2),4),";  ",
        "V2: ",myRound(midPoint(V2),4),";  ",
        "A: " ,myRound(midPoint(A),4),sep="")
    str[2]=paste("called ",CALLED,"; depth ",depth,"; max ",MAXDEPTH,"; ")
    str[3]=paste("strong ",STRONG,"; weak ",WEAK,"  ")
    str[4]=paste("passed ",duration2string(passed),"; total h ",total,"; done ",done)
    cat("\r")
    for (i in 1:4){
      cat(str[i])
      w=getOption("width")-nchar(str[i])+4
      w=max(w,0)
      cat(rep(" ",w),sep="")
    }
  }
  
  
  
  T_ID[ID]<<-ID
  
  T_T1_min[ID]<<-T1[1]
  T_T1_max[ID]<<-T1[2]
  T_V1_min[ID]<<-V1[1]
  T_V1_max[ID]<<-V1[2]
  T_T2_min[ID]<<-T2[1]
  T_T2_max[ID]<<-T2[2]
  T_V2_min[ID]<<-V2[1]
  T_V2_max[ID]<<-V2[2]
  T_A_min[ID] <<-A[1]
  T_A_max[ID] <<-A[2]
  

  #1) try to apply linear Theorem
  
  masterAnalyserWrapper(Points=Points,T1=T1,V1=V1,T2=T2,V2=V2,A=A)
  
  for (i in 1:min(nrow(master),4)){
    if (weakTheorem(Points=Points,T1=T1,V1=V1,T2=T2,V2=V2,A=A,NNN=i)){
      T_nodetype[ID]<<-1
      T_wx[ID]<<-global[1]
      T_wy[ID]<<-global[2]
      T_S_index[ID]<<-global[3]
      
      
      done<<-done+percentage
      WEAK<<-WEAK+1
      return(TRUE)
    }
  }
  
  
  #2) try to apply strong Theorem
  if (strongTheorem(Points=Points,T1=T1,V1=V1,T2=T2,V2=V2,A=A)){
    T_nodetype[ID]<<-2
    
    T_P1_index[ID]<<-global[1]
    T_P2_index[ID]<<-global[2]
    T_P3_index[ID]<<-global[3]
    T_Q1_index[ID]<<-global[4]
    T_Q2_index[ID]<<-global[5]
    T_Q3_index[ID]<<-global[6]
    T_s_p[ID]<<-global[7]  ## will eventually be unused
    T_s_q[ID]<<-global[8]
    T_r[ID]<<-global[9]
    
    done<<-done+percentage
    STRONG<<-STRONG+1
    return(TRUE)
  }
  
  
  #3) We need to go deeper
  
  if (depth<=4){
    t1_splits=1
    v1_splits=1
    t2_splits=1
    v2_splits=1
    a_splits =1
    if (depth==0){T_split[ID]<<-1;t1_splits=4}
    if (depth==1){T_split[ID]<<-2;v1_splits=30}
    if (depth==2){T_split[ID]<<-3;t2_splits=4}
    if (depth==3){T_split[ID]<<-4;v2_splits=15}
    if (depth==4){T_split[ID]<<-5;a_splits =30}
  } else {
    t1_splits=2
    v1_splits=2
    t2_splits=2
    v2_splits=2
    a_splits =2
    T_split[ID]<<-6
  }

  
  T_nodetype[ID]<<-3
  T_nrChildren[ID]<<-t1_splits*v1_splits*t2_splits*v2_splits*a_splits
  
  IDfirstChild=nextUnused
  T_IDfirstChild[ID]<<-nextUnused
  nextUnused<<-nextUnused+T_nrChildren[ID]
  
  childCounter=0
  for (t1_index in 1:t1_splits){
    for (v1_index in 1:v1_splits){
      for (t2_index in 1:t2_splits){
        for (v2_index in 1:v2_splits){
          for (a_index in 1:a_splits){
            res=RupertDisprover(Points=Points,
                                T1=splitInterval(T1,t1_index,t1_splits),
                                V1=splitInterval(V1,v1_index,v1_splits),
                                T2=splitInterval(T2,t2_index,t2_splits),
                                V2=splitInterval(V2,v2_index,v2_splits),
                                A =splitInterval(A , a_index, a_splits),
                                depth=depth+1,
                                ID=IDfirstChild+childCounter)
            
            childCounter=childCounter+1
            
            if (res==FALSE){return(FALSE)}
            
          }
        }
      }
    }
  }
  
  return(TRUE)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 7. Loading the solid ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~


## this is the new sexy polyhedron:
C1=c(152024884,0,210152163 )/259375205
C2=c(0.6632738028,0.6106948881,0.3980949609)
C3=c(0.8193990033,0.5298215096,0.1230614493)



Points=matrix(nrow=90,ncol=3)
for (i in 0:2){
  for (k in 0:14){
    for (l in 0:1){
      if (i==0){Points[k+15*i+45*l+1,]=(-1)^l*Rz(2*pi*k/15)%*%C1}
      if (i==1){Points[k+15*i+45*l+1,]=(-1)^l*Rz(2*pi*k/15)%*%C2}
      if (i==2){Points[k+15*i+45*l+1,]=(-1)^l*Rz(2*pi*k/15)%*%C3}
    }
  }
}
rm(C1,C2,C3)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 8. Generating the solution tree ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


n=19*10^6 ## that much space will be allocated
T_ID=rep(NA,n)
T_nodetype=rep(NA,n)
T_nrChildren=rep(NA,n)
T_IDfirstChild=rep(NA,n)
T_split=rep(NA,n)

T_T1_min=rep(NA,n)
T_T1_max=rep(NA,n)
T_V1_min=rep(NA,n)
T_V1_max=rep(NA,n)
T_T2_min=rep(NA,n)
T_T2_max=rep(NA,n)
T_V2_min=rep(NA,n)
T_V2_max=rep(NA,n)
T_A_min =rep(NA,n)
T_A_max =rep(NA,n)

T_wx =rep(NA,n)
T_wy =rep(NA,n)
T_S_index =rep(NA,n)

T_P1_index =rep(NA,n)
T_P2_index =rep(NA,n)
T_P3_index =rep(NA,n)
T_Q1_index =rep(NA,n)
T_Q2_index =rep(NA,n)
T_Q3_index =rep(NA,n)

T_r =rep(NA,n)
T_s_p =rep(NA,n) ##will be always 1
T_s_q =rep(NA,n)
nextUnused=2



RupertDisprover(Points,
                T1=c(0,0.42),V1=c(0,3.15),
                T2=c(0,0.42),V2=c(0,1.58),
                A=c(-1.58,1.58))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 9. Processing and storing the data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df=data.frame(T_ID,T_nodetype,T_nrChildren, ## meta stuff
              T_IDfirstChild,T_split,
              T_T1_min,T_T1_max,T_V1_min,T_V1_max, ## intervals
              T_T2_min,T_T2_max,T_V2_min,T_V2_max,
              T_A_min ,T_A_max,
              T_wx,T_wy,T_S_index, ## weak shit
              T_P1_index,T_P2_index,T_P3_index,## strong shit
              T_Q1_index,T_Q2_index,T_Q3_index,
              T_r,T_s_p,T_s_q)
df=as.data.table(df)

df=df[1:(nextUnused-1)]

colnames(df)=substr(colnames(df),3,100)


## Step 1) Storing the interval endpoints as integers


factorIntervals=2^13*3*5^4 ## = 15360000 seems to be optimal

is.int(df[,T1_min]*factorIntervals)
is.int(df[,T1_max]*factorIntervals)
is.int(df[,V1_min]*factorIntervals)
is.int(df[,V1_max]*factorIntervals)
is.int(df[,T2_min]*factorIntervals)
is.int(df[,T2_max]*factorIntervals)
is.int(df[,V2_min]*factorIntervals)
is.int(df[,V2_max]*factorIntervals)
is.int(df[,A_min ]*factorIntervals)
is.int(df[,A_max ]*factorIntervals)


## using data.table syntax
df[,T1_min:=round(T1_min*factorIntervals)]
df[,T1_max:=round(T1_max*factorIntervals)]
df[,V1_min:=round(V1_min*factorIntervals)]
df[,V1_max:=round(V1_max*factorIntervals)]
df[,T2_min:=round(T2_min*factorIntervals)]
df[,T2_max:=round(T2_max*factorIntervals)]
df[,V2_min:=round(V2_min*factorIntervals)]
df[,V2_max:=round(V2_max*factorIntervals)]
df[,A_min :=round(A_min *factorIntervals)]
df[,A_max :=round(A_max *factorIntervals)]




## Step 2) Making w in the weak Theorem rational

## we don't want to access the df unnecessarily 
wx=df$wx
wy=df$wy

## here we will store the results:
wx_nominator=rep(NA,nrow(df))
wy_nominator=rep(NA,nrow(df))
w_denominator=rep(NA,nrow(df))

for (i in which(df$nodetype==1)){
  if (i%%1000==0){cat(i,"\r")}
  vec=c(wx[i],wy[i])
  vecRational=approx_on_circle(vec)
  
  wx_nominator[i]=vecRational[1]
  wy_nominator[i]=vecRational[2]
  w_denominator[i]=vecRational[3]
} ## takes 15-30 minutes


df$wx_nominator=wx_nominator
df$wy_nominator=wy_nominator
df$w_denominator=w_denominator

## Step 3) Removing unnecessary stuff:

df$wx=NULL
df$wy=NULL
table(df$s_p) ## is always 1
df$s_p=NULL ## removing s_p


## Step 4) Saving everything:

## 3 ways to save the same thing:

## 1)
saveRDS(df,"sexypoly-solutionTree_try22012025.rds")

## 2) FUCK FWRITE
data.table::fwrite(df,"sexypoly-solutionTree_try22012025.csv",row.names = F,sep = ";")


