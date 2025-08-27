from sage.all import *

def sinQ(x):
    assert isinstance(x, Rational) or isinstance(x, Integer), f"Give me rational inputs, please! You gave me {type(x)}"
    assert abs(x)<=4, f"x should be between -4 and 4 and you gave me {x}"
    
    x2 = x**2
    
    return x*(1-x2/(2*3)*
             (1-x2/(4*5)*
             (1-x2/(6*7)*
             (1-x2/(8*9)*
             (1-x2/(10*11)*
             (1-x2/(12*13)*
             (1-x2/(14*15)*
             (1-x2/(16*17)*
             (1-x2/(18*19)*
             (1-x2/(20*21)*
             (1-x2/(22*23)*
             (1-x2/(24*25)))))))))))))

def cosQ(x):
    assert isinstance(x, Rational) or isinstance(x, Integer), f"Give me rational inputs, please! You gave me {type(x)}"
    assert abs(x)<=4, f"x should be between -4 and 4 and you gave me {x}"
    
    x2 = x**2
  
    return (1-x2/(1*2)*
           (1-x2/(3*4)*
           (1-x2/(5*6)*
           (1-x2/(7*8)*
           (1-x2/(9*10)*
           (1-x2/(11*12)*
           (1-x2/(13*14)*
           (1-x2/(15*16)*
           (1-x2/(17*18)*
           (1-x2/(19*20)*
           (1-x2/(21*22)*
           (1-x2/(23*24)))))))))))))

def ScalarProduct(vector1,vector2):
    ## returns the scalar product of two (vertical) vectors, aka n times 1 matrices
    
    assert vector1.ncols()==1, f"ScalarProduct wants integers or matrix as input"
    assert vector2.ncols()==1, f"ScalarProduct wants integers or matrix as input"
    assert vector1.nrows()==vector2.nrows(), f"Matrices must have the same number of rows"
    
    assert vector1.base_ring()==ZZ or vector1.base_ring()==QQ, f"ScalarProduct wants integers or matrix as input"
    assert vector2.base_ring()==ZZ or vector2.base_ring()==QQ, f"ScalarProduct wants integers or matrix as input"

    return (vector1.T*vector2)[0,0]

def NormSquared(vector):
    ## returns the squared norm of a (vertical) vector, an nx1 matrix
    return ScalarProduct(vector,vector)

def abs_matrix(M): 
    ## every entry is replaced by its absolute value
    nrows = M.nrows()
    ncols = M.ncols()
    new_matrix = matrix(M.base_ring(), nrows, ncols) #Create an empty matrix with the same base ring
    for i in range(nrows):
        for j in range(ncols):
            new_matrix[i,j] = abs(M[i,j])
    return new_matrix

def R(alpha):
    assert isinstance(alpha, Rational) or isinstance(alpha, Integer), f"Give me rational inputs, please! You gave me {type(alpha)}"

    A=matrix(QQ,2,2)
    A[0,:]=matrix([cosQ(alpha),-sinQ(alpha)])
    A[1,:]=matrix([sinQ(alpha),cosQ(alpha)])
    return A

def Rz(alpha):
    return matrix([[cos(alpha),-sin(alpha),0],[sin(alpha),cos(alpha),0],[0,0,1]])

def X(theta,phi):
    assert isinstance(theta, Rational) or isinstance(theta, Integer), f"Give me rational inputs, please! You gave me {type(theta)}"
    assert isinstance(phi  , Rational) or isinstance(phi  , Integer), f"Give me rational inputs, please! You gave me {type(phi  )}"
    
    A=matrix(QQ,3,1)
    A[:,0]=matrix([[cosQ(theta)*sinQ(phi)],[sinQ(theta)*sinQ(phi)],[cosQ(phi)]])

    return A 

def M(theta,phi):
    assert isinstance(theta, Rational) or isinstance(theta, Integer), f"Give me rational inputs, please! You gave me {type(theta)}"
    assert isinstance(phi  , Rational) or isinstance(phi  , Integer), f"Give me rational inputs, please! You gave me {type(phi  )}"
    
    A=matrix(QQ,2,3)
    A[0,:]=matrix([-sinQ(theta),cosQ(theta),0])
    A[1,:]=matrix([-cosQ(theta)*cosQ(phi),-sinQ(theta)*cosQ(phi),sinQ(phi)])

    return A 

def M_theta_prime(theta,phi):
    assert isinstance(theta, Rational) or isinstance(theta, Integer), f"Give me rational inputs, please! You gave me {type(theta)}"
    assert isinstance(phi  , Rational) or isinstance(phi  , Integer), f"Give me rational inputs, please! You gave me {type(phi  )}"
    
    A=matrix(QQ,2,3)
    A[0,:]=matrix([-cosQ(theta),-sinQ(theta),0])
    A[1,:]=matrix([sinQ(theta)*cosQ(phi),-cosQ(theta)*cosQ(phi),0])

    return A

def M_phi_prime(theta,phi):
    assert isinstance(theta, Rational) or isinstance(theta, Integer), f"Give me rational inputs, please! You gave me {type(theta)}"
    assert isinstance(phi  , Rational) or isinstance(phi  , Integer), f"Give me rational inputs, please! You gave me {type(phi  )}"
    
    A=matrix(QQ,2,3)
    A[0,:]=matrix([0,0,0])
    A[1,:]=matrix([cosQ(theta)*sinQ(phi),sinQ(theta)*sinQ(phi),cosQ(phi)])

    return A

def R_alpha_prime(alpha):
    assert isinstance(alpha, Rational) or isinstance(alpha, Integer), f"Give me rational inputs, please! You gave me {type(alpha)}"

    A=matrix(QQ,2,2)
    A[0,:]=matrix([-sinQ(alpha),-cosQ(alpha)])
    A[1,:]=matrix([cosQ(alpha),-sinQ(alpha)])

    return A
    
def sqrt_lower(x):
    ## returns a non-negative rational value, that is smaller (or equal) to the square root of the input
    
    assert isinstance(x, Rational) or isinstance(x, Integer), f"Give me rational inputs, please! You gave me {type(x)}"
    assert x>=0, f"I really tried to find the square root of {x}, but I couldn't"

    if x==0:
        return 0
    
    # first we find the unique integer a such that x*100**a in [10**20,10**22)                  

    ## initial guess for a:
    a=-floor(log(float(x))/log(100.))+10
  
    while x*100**a < 10**20 or x*100**a >= 10**22: ## repeat until x*100**a in [10**20,10**22)
        if x*100**a < 10**20:
            a+=1
        if x*100**a >= 10**22:
            a-=1

    ## we now have x*100**a in [10**20,10**22)

    ## we find the unique(!) positive integer b such that b**2<= x*100**a < (b+1)**2 
      
    ## initial guess for b:
    b=floor(sqrt(float(x*100**a)))
  
    while b**2> x*100**a or x*100**a >= (b+1)**2:
        if b**2 > x*100**a:
            b-=1
        if (b+1)**2 <= x*100**a:
            b+=1
            
    assert b>0, "b should be positive"
            
    result=b/10**a
    
    ## a few asserts to ensure that everything works properly: 
    assert result>=0, f"The approximation should be non-negative"
    assert result**2<=x, f"We got {result} as a lower approximation of the sqrt of {x}"
    assert isinstance(result, Rational) or isinstance(result, Integer), f"Something went wrong!"
    
    return result

def sqrt_upper(x):
    ## returns a non-negative rational value, that is greater (or equal) to the square root of the input

    assert isinstance(x, Rational) or isinstance(x, Integer), f"Give me rational inputs, please! You gave me {type(x)}"
    assert x>=0, f"I really tried to find the square root of {x}, but I couldn't"

    if x==0:
        return 0
    
    result=x/sqrt_lower(x)

    assert result>=0, f"The approximation should be non-negative"
    assert result**2>=x, f"We got {result} as an upper approximation of the sqrt of {x}"
    assert isinstance(result, Rational) or isinstance(result, Integer), f"Something went wrong!"
    
    return result

def R_exact(alpha): # Needed for the Ruperthedron
    A=matrix(SR,2,2)
    A[0,:]=matrix([cos(alpha),-sin(alpha)])
    A[1,:]=matrix([sin(alpha),cos(alpha)])
    return A

def M_exact(theta,phi): # Needed for the Ruperthedron
    A=matrix(SR,2,3) 
    A[0,:]=matrix([-sin(theta),cos(theta),0])
    A[1,:]=matrix([-cos(theta)*cos(phi),-sin(theta)*cos(phi),sin(phi)])
    return A

def ScalarProduct_exact(vector1,vector2): # Needed for the Ruperthedron
    ## returns the scalar product of two (vertical) vectors, aka n times 1 matrices
    
    assert vector1.ncols()==1, f"ScalarProduct wants integers or matrix as input"
    assert vector2.ncols()==1, f"ScalarProduct wants integers or matrix as input"
    assert vector1.nrows()==vector2.nrows(), f"Matrices must have the same number of rows"
    
    return (vector1.T*vector2)[0,0]