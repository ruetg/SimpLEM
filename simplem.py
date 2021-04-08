from numba import jit
from numba.experimental import jitclass
from numba import int64,float64
import numpy
import matplotlib.pyplot as plt
import math
import timeit
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter



spec2= [
        ('__nn', int64),               
        ('__numel',int64),
        ('__u', int64),
        ('__uu', int64),
        ('__ul', int64),
        ('__left', int64[:]),
        ('__z', float64[:])
    ]
@jitclass(spec2)
    
class pq:
        def __init__(self,z):
            
            self.__nn = numpy.int64(len(z))
            self.__numel = numpy.int64(0)
            self.__u = numpy.int64(0)
            self.__uu = numpy.int64(0)
            self.__ul = numpy.int64(0)
            self.__left = numpy.full(len(z)+1,0)
            self.__z = numpy.concatenate((numpy.zeros(1).ravel(),z.ravel()))
        def top(self):
            return self.__left[1]-1
        def get(self):
            return self.__z[self.__left]
        def pop(self):
            self.__uu = self.__left[1]
            self.__left[1] = self.__left[self.__numel]
            self.__left[self.__numel] = 0
            self.__u = 2
            self.__ul = numpy.int(self.__u/2)
            while (self.__u < self.__numel - 2):
                if (self.__z[self.__left[self.__u]] < self.__z[self.__left[self.__u + 1]]):

                    if (self.__z[self.__left[self.__ul] ]> self.__z[self.__left[self.__u]]):
                        t = self.__left[self.__ul]
                        self.__left[self.__ul] = self.__left[self.__u]
                        self.__left[self.__u] = t

                        self.__ul = self.__u
                        self.__u *= 2
                    else:
                        break;
                    
                elif(self.__z[self.__left[self.__ul]] > self.__z[self.__left[self.__u + 1]]):
                
                    t = self.__left[self.__ul]
                    self.__left[self.__ul] = self.__left[self.__u + 1]
                    self.__left[self.__u + 1] = t
                    self.__u= 2*(self.__u+1);
                    self.__ul=numpy.int(self.__u / 2);
                
                else:
                
                    break;

            self.__numel -= 1;
            return self;

        def push(self, i):
            i+=1
            self.__numel+=1


            self.__u = self.__numel;
            self.__ul =  numpy.int(self.__u / 2);

            self.__left[self.__u] = i;
            while (self.__ul > 0):
                if (self.__z[self.__left[self.__ul]] > self.__z[self.__left[self.__u]]):

                    t = self.__left[self.__ul]
                    self.__left[self.__ul] = self.__left[self.__u]
                    self.__left[self.__u] = t

                else:
                    break;
                self.__u=numpy.int(self.__u/2)
                self.__ul=numpy.int(self.__u/2)
            return self




if 1:
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
        ('U',float64),
#         ('BCX',int64[:,:]),
#         ('BC',int64[:])

    ]
    

@jitclass(spec)
class fs:
    def __init__(self):
        #Modify parameters here
        self.m=.5 #Drainage area exponent
        self.dx=1000.0 # grid spacing (m)
        self.dy=1000.0 
        self.t=100e6 #total time (yr)
        self.dt=1e6 # Time step
        self.nx=499 #Number of x grid points
        self.ny=500 #Number of y grid points
        self.A=numpy.ones((self.ny,self.nx),dtype=numpy.float64)
        self.Z=numpy.random.rand(self.ny,self.nx) *10
        self.k=1e-6 #Erodibility 
        self.U=.0001 #Uplift rate
#         self.BCX = numpy.full(numpy.shape(self.Z),0)
#         self.BCX[:,0] = 1
#         self.BCX[:,-1] = 1
#         self.BCX[0,:] = 1
#         self.BCX[-1,:] = 1
#         self.BC = numpy.where(self.BCX == 1)[0]
        self.s=numpy.zeros((self.ny,self.nx),dtype=numpy.int64)
        self.I=numpy.zeros(self.ny*self.nx,dtype=numpy.int64)
    def sinkfill(self):
            pittop = 0
            c=int(0)
            nn = self.nx*self.ny
            p = int(0)
            closed = numpy.full(nn,False)
            pittop = numpy.zeros(nn)
            pit = numpy.zeros(nn)
            idx = [1,-1,self.ny,-self.ny,-self.ny+1,-self.ny-1,self.ny+1,self.ny-1]
            open = pq(self.Z.transpose().flatten());
            
           # for i in range(len(self.BC)):
                #open = open.push(self.BC[i]);
                #closed[self.BC[i]]=True;

            for i in range(0,self.ny):
                if (closed[i]==False):
                
                    closed[i]=True;

                    open = open.push(i);

                    c+=1;
                    
            for  i in range(self.ny*self.nx-self.ny, self.ny*self.nx):
                if (closed[i]==False):
                    closed[i]=True;
                    open = open.push(i);
                    c+=1;
            for  i in range(self.ny-1,self.nx*self.ny,self.ny):
                if (closed[i]==False):
                    closed[i]=True;
                    open = open.push(i);
                    c+=1;
            for  i in range(self.ny,nn,self.ny):
                if (closed[i]==True):
                    closed[i]=True

                    open = open.push(i)

                    c+=1;
            s=int(0)
            si = int(0);
            ij  = int(0)
            ii = int(0)
            jj = int(0)
            ci = int(0)
            pittop=int(-9999)
            while ((c>0) or (p>0)):
                if ((p>0) and (c>0) and (pit[p-1] == -9999)):
                    s = open.top()
                    open = open.pop();
                    c-=1
                    pittop=-9999
                elif (p>0):
                    s=int(pit[p-1])
                    pit[p-1]=-9999;
                    p-=1
                    if (pittop==-9999):
                        si,sj = self.lind(s,self.ny)
                        pittop=self.Z[si,sj]
                else:
                    s=int(open.top());
                    open=open.pop();
                    c-=1
                    pittop=-9999;
                    
                for i in range(8):

                    ij = idx[i]+s;
                    si,sj = self.lind(s,self.ny)
                    ii,jj = self.lind(ij,self.ny)
                    if ((ii>=0) and (jj>=0) and (ii<self.ny) and (jj<self.nx)):
                        if (closed[ij]==False):
                            closed[ij]=True;


                            if (self.Z[ii,jj]<=self.Z[si,sj]):

                                self.Z[ii,jj]=self.Z[si,sj]+1e-7;

                                pit[p]=ij
                                p+=1
                            else:
                                open = open.push(ij);
                                c+=1

    def lind(self,xy,n):#Compute linear index from 2 point
        x=math.floor(xy/n)
        y=xy%n
        return y,x
    def set_z(self,Z):
        ny,nx = np.shape(Z)
        print('here')
        self.Z = Z
        self.A=numpy.ones((self.ny,self.nx),dtype=numpy.float64)
      
        self.s=numpy.zeros((self.ny,self.nx),dtype=numpy.int64)
        self.I=numpy.zeros(self.ny*self.nx,dtype=numpy.int64)
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
                            
if __name__ == '__main__':
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

                        
    



