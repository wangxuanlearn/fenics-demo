from fenics import*
from mshr import *
import matplotlib.pyplot as plt
domain=Circle(Point(0,0),1)
mesh=generate_mesh(domain,64)
V=FunctionSpace(mesh,'P',1)
u_D=Expression('0',degree=2)
def boundary(x,on_boundary):
    print(on_boundary)
    print(type(on_boundary))
    return on_boundary
#print(u_D)
bc=DirichletBC(V,u_D,boundary)

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
plot(p,title='Load')           
vtkfile_w=File('poisson_membrane/deflection.pvd')
vtkfile_w<<w

vtkfile_p=File('poisson_membrane/load.pvd')
vtkfile_p<<p