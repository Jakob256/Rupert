Rz <- function(a){
  return(matrix(c(cos(a),sin(a),0,-sin(a),cos(a),0,0,0,1),3,3))
}

## this is the new sexy polyhedron:
C1=c(152024884,0,210152163 )/259375205
C2=c(0.6632738028,0.6106948881,0.3980949609)
C3=c(0.8193990033,0.5298215096,0.1230614493)

C1=c(3939,0,4340)/5861
C2=c(7855,4178,4484)/10^4
C3=c(9526,2057,1102)/10^4

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

ds = c()
ds = c(ds, sqrt(sum((Points[1,] - Points[2,])^2)))
ds = c(ds, sqrt(sum((Points[1+15,] - Points[2+15,])^2)))
ds = c(ds, sqrt(sum((Points[1+30,] - Points[2+30,])^2)))

ds

dtemp = c()
for(j in 1:15){
  for(i in 16:30){
    dtemp = c(dtemp,sqrt(sum((Points[i,] - Points[j,])^2)))
  }
}
ds = c(ds,unique(round(sort(dtemp),5))[1:2])
ds

dtemp = c()
for(j in 16:30){
  for(i in 31:45){
    dtemp = c(dtemp,sqrt(sum((Points[i,] - Points[j,])^2)))
  }
}
ds = c(ds,unique(round(sort(dtemp),5))[1:2])

dtemp = c()
for(j in 31:45){
  for(i in 76:90){
    dtemp = c(dtemp,sqrt(sum((Points[i,] - Points[j,])^2)))
  }
}

ds = c(ds,unique(round(sort(dtemp),5))[1])



library(rgl)

vertices = t(Points)

indices = matrix(nrow = 3, ncol = 0)

ds

rds = round(ds,3)

for(j in 1:90){
  for(i in 1:90){
    for(k in 1:90){
      for(d in ds){
        dd1 = round(sqrt(sum((Points[i,] - Points[j,])^2)),3)
        dd2 = round(sqrt(sum((Points[j,] - Points[k,])^2)),3)
        dd3 = round(sqrt(sum((Points[k,] - Points[i,])^2)),3)
        if(any(dd1 == rds) && any(dd2 == rds) && any(dd3 == rds)){
          indices = cbind(indices, c(i,j,k))
        }
      }
    }
  }
}

indices_bak = indices

N1 = dim(indices_bak[, !duplicated(t(indices_bak))])[2]

for(i in 2:14){
  indices_bak = cbind(indices_bak, c(1,i,i+1))
#  indices_bak = cbind(indices_bak, c(2,2*i,2*(i+1)))
}

indices_unique <- indices_bak[, !duplicated(t(indices_bak))]

bg3d(color = "white")  # Change background to white

# Clear previous lights and add new light with soft shadows
clear3d(type = "lights")
light3d(theta = 30, phi = 30, specular = "grey70", diffuse = "grey80", ambient = "grey40")

# Adjust the material properties for a less shiny and more diffused surface
#material3d(
#  color = "red",     # Softer grey color
#  ambient = "grey50",   # Ambient light effect
#  specular = "grey40",  # Soft reflection
#  shininess = 10,       # Reduce shininess for a matte look
#  emission = "grey70"    # Remove glow
#)

# Create a 3D mesh object from vertices and indices
solid <- tmesh3d(
  vertices = vertices,  # Matrix of vertices
  indices = indices_unique,    # Matrix of indices
  homogeneous = FALSE   # Not using homogeneous coordinates
)

indices_unique

# Render the solid object with transparency and reduced brightness
shade3d(solid, color = "grey80", alpha = 1)  # Slightly darker and less transparent

N = dim(indices_unique)[2]

for(i in 1:N1){
  lines3d(t(vertices[,indices_unique[,i]]), color = "black", lwd = 2)
}

clear3d(type = "lights")
#light3d(theta = 60, phi = 30, viewpoint.rel = T, specular = "grey100", diffuse = "grey100", ambient = "grey100")
#light3d(theta = -10, phi = -10, viewpoint.rel = T, specular = "grey60", diffuse = "grey60", ambient = "grey60")

light3d(theta = 60, phi = 30, specular = "grey70", diffuse = "grey80", ambient = "grey40")
light3d(theta = -10, phi = -10, specular = "grey70", diffuse = "grey80", ambient = "grey40")

rgl.snapshot("Nolocalhedron.png", fmt = "png")

#quads3d(x, y, z, color = "lightblue", alpha = 0.4)

# Create a mesh3d plane at z = -1
plane_mesh <- tmesh3d(
  vertices = rbind(
    x, y, z  # z all = -1
  ),
  indices = matrix(c(1,2,3, 1,3,4), 3, 2),
  homogeneous = FALSE
)
#shade3d(plane_mesh, color="lightblue", alpha=0.4)

shadow3d(obj = plane_mesh, mesh = solid,
         up = c(0, 0, 1),
         col = "black",
         alpha = 1,
         minVertices = 10000)

