from fenics import*
#from mshr import *
#import mshr
import mshr as mr
import numpy as np
import matplotlib.pyplot as plt
domain=mr.Circle(Point(0,0),1)
mesh=mr.generate_mesh(domain,64)
V=FunctionSpace(mesh,'P',1)
w_D=Expression('0',degree=2)
def boundary(x,on_boundary):
    return on_boundary
bc=DirichletBC(V,w_D,boundary)
beta=8
R0=0.6
p=Expression('4*exp(-pow(beta,2)*(pow(x[0],2)+pow(x[1]-R0,2)))',
              degree=1,beta=beta,R0=R0)
w=TrialFunction(V)
v=TestFunction(V)
a=dot(grad(w),grad(v))*dx
L=p*v*dx  
w=Function(V) 
solve(a==L,w,bc) 
p=interpolate(p,V) 
plot(w,title='Deflection') 
plt.show()          
vtkfile_w=File('poisson_membrane/deflection.pvd')
vtkfile_w << w

plot(p,title='Load') 
vtkfile_p=File('poisson_membrane/load.pvd')
vtkfile_p << p
plt.show()
tol=0.001 #avoid hitting points outside the domain
y=np.linspace(-1+tol,1-tol,101)
points=[(0,y_)for y_ in y]#2D points
w_line=np.array([w(point) for point in points])
p_line=np.array([p(point) for point in points])
# print(50*w_line)
# print(p_line)
plt.plot(y,50*w_line,'k',linewidth=2)
plt.plot(y,p_line,'b--',linewidth=2)
plt.grid(True)
plt.xlabel('$y$')
# plt.xlim((-1,1))
plt.legend(['Deflection($\\times 50$)','Load'],loc='upper left')
plt.savefig('poisson_membrane/curves.pdf')
plt.savefig('poisson_membrane/curves.png')
plt.show()