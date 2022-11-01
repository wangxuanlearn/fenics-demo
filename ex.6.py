from fenics import *
import numpy as np
import sympy as sym
#Define nonlinear coefficient
def q(u):
    "Return nonlinear coefficient"
    return 1 +u**2

#Create mesh and define function space
mesh=UnitSquareMesh(8,8)
V=FunctionSpace(mesh,'P',1)

#Use SymPy to compute f from the manufactured solution u
x, y = sym.symbols('x[0],x[1]')
u = 1+x+2*y
f = -sym.diff(q(u)*sym.diff(u,x),x)-sym.diff(q(u)*sym.diff(u,y),y)
f = sym.simplify(f)
u_code = sym.printing.ccode(u)  
f_code = sym.printing.ccode(f)
print('u=', u_code)
print('f=', f_code)
#T = 2.0
#num_steps = 50 # number of time steps
#dt=T / num_steps #time step size
#alpha=3
#beta=1.2

#nx = ny = 30
#mesh=RectangleMesh(Point(-2,-2),Point(2,2),nx,ny)

#Define boundary condition
u_D = Expression(u_code , degree=1)

def boundary(x , on_boundary):
    return on_boundary
bc = DirichletBC(V, u_D, boundary)

#Define initial value
#u_0=Expression('exp(-a*pow(x[0],2)-a*pow(x[1],2))',degree=2,a=5)
#u_n=interpolate(u_0,V)

#Define variational problem
u = Function(V)
v = TestFunction(V)
f = Expression(f_code, degree=1)
#u_n=interpolate(u_0,V)


#Define variational problem
F=q(u)*dot(grad(u),grad (v))*dx-f*v*dx
solve(F==0, u, bc)
#a,L=lhs(F),rhs(F)

#Plot solution
plot(u,title='Velocity')


 #Compute error at vertices
u_e=interpolate(u_D,V)
error_max=np.abs(u_e.vector()-u.vector()).max()
print('error_max=', error_max)


#Create file for saving solution
#ufile = File("heat-guassian/solution.pvd")

#Time-stepping

#t=0
#for n in range(num_steps):
    #Update current time
    #t +=dt
    
    #Solve variational problem
    #solve(a==L, u, bc)

    
    #ufile  << (u,t)
   
    #Update previous solution
    #u_n.assign(u)