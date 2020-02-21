from numba import jit,jitclass
from numba import int64,float64
import numpy
import matplotlib.pyplot as plt
import math
import timeit
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
spec = [
    ('m', float64),               
    ('dx',float64),
    ('dy',float64),
    ('t',float64),
    ('dt',float64),
    ('nx',int64),
    ('ny',int64),
    ('A', float64[:,:]),
    ('Z',float64[:,:]),
    ('k',float64),
    ('m',float64),
    ('s',int64[:,:]),
    ('I',int64[:]),         
    ('U',float64)

]

@jitclass(spec)
class fs:
    def __init__(self):
        #Modify parameters here
        self.m=.5 #Drainage area exponent
        self.dx=1000.0 # grid spacing (m)
        self.dy=1000.0 
        self.t=100e6 #total time (yr)
        self.dt=1e5 # Time step
        self.nx=500 #Number of x grid points
        self.ny=500 #Number of y grid points
        self.A=numpy.ones((self.ny,self.nx),dtype=numpy.float64)
        self.Z=numpy.random.rand(self.ny,self.nx) 
        self.k=1e-6 #Erodibility 
        self.U=.0001 #Uplift rate

        self.s=numpy.zeros((self.nx,self.ny),dtype=numpy.int64)
        self.I=numpy.zeros(self.ny*self.nx,dtype=numpy.int64)
    def lind(self,xy,n):#Compute linear index from 2 point
        x=math.floor(xy/n)
        y=xy%n
        return y,x
    def slp(self):#Calculate slope and steepest descent
        ij=0  
        for i in range(0,self.ny):
            for j in range(0,self.nx):
                ij=j*self.ny+i
                mxi=0
                self.s[i,j]=ij
                if (i>0 and i<self.ny and j>0 and j<self.nx-1 and i<self.ny-1):
                    for i1 in range(-1,2):
                        for j1 in range(-1,2):
                            if (self.Z[i,j]-self.Z[i+i1,j+j1])/numpy.sqrt((float(i1)**2)+float(j1)**2+1e-10)>mxi:
                                ij2=(j+j1)*self.ny+i1+i
                                mxi=(self.Z[i,j]-self.Z[i+i1,j+j1])/numpy.sqrt((float(i1)**2)+float(j1)**2)
                                self.s[i,j]=ij2
    def stack(self): #Generate the stack
        c=0
        k=0
        for i in range(0,self.ny):
            for j in range(0,self.nx):

                ij=j*self.ny+i
                i2=i
                j2=j
                if self.s[i,j]==ij:
                    
                    self.I[c]=ij
                    c+=1
                    
                    while c>k and c<self.ny*self.nx-1:
                        for i1 in range(-1,2):
                            for j1 in range(-1,2):
                                if j2+j1>0 and i2+i1>0 and j2+j1<self.nx-1 and i2+i1<self.ny-1:
                                    ij2=(j2+j1)*self.ny+i2+i1
                                    if ij!=ij2 and self.s[i2+i1,j2+j1]==ij:
                                        self.I[c]=ij2
                                        c+=1
                                        
                        
                        k=k+1
                        ij=self.I[k]
                        i2,j2=self.lind(ij,self.ny)


    def acc(self):#Calculate drainage area
        self.A[:,:]=1.0
        for ij in range(len(self.I)-1,0,-1):
            i,j=self.lind(self.I[ij],self.ny)
            i2,j2=self.lind(self.s[i,j],self.ny)
            self.A[i2,j2]+=self.A[i,j]
    def erode(self):#Erode using fastscape method
        dA=(self.dx*self.dy)**self.m
        f=self.dt/self.dx*self.k*dA
        for ij in range(0,len(self.I)):
            i,j=self.lind(self.I[ij],self.ny)
            i2,j2=self.lind(self.s[i,j],self.ny)
            self.Z[i,j]=1.0/(1.0+f*self.A[i,j]**self.m)*(self.A[i,j]**self.m*f*self.Z[i2,j2]+self.Z[i,j])
        self.Z+=self.U*self.dt
        self.Z[:,-1]=0
        self.Z[:,0]=0
        self.Z[-1,:]=0
        self.Z[0,:]=0
    #def diffusion(self):
     #   for i in range(2,self.nx-1):
      #      for j in range(2,self.ny-1):
       #         self.Z(i,j)=k*(self.Z(i,j+1)+self.Z(i,j)*2
                            

A=fs()
X,Y=numpy.meshgrid(range(0,A.nx),range(0,A.ny))
fig=plt.figure()

for t in range(0,int(A.t/A.dt)): #main loop
    start=timeit.default_timer()

    A.slp()
    A.stack()
    A.acc()
    A.erode()
    end=timeit.default_timer()
    print(end-start)

    a=plt.imshow(A.Z)
    plt.colorbar(a)

    plt.pause(.05)
    plt.clf()

    
print(numpy.where(A.I < 1))

                        






