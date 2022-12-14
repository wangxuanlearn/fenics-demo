
from __future__ import print_function
from fenics import *
#from mshr import *
import numpy as np
T = 10.0           # final time
num_steps=500
dt=T/num_steps
mu=1.0
rho=1.0

#Create mesh and define function space
mesh=UnitSquareMesh(16,16)
V=VectorFunctionSpace(mesh,'P',2)
Q=FunctionSpace(mesh,'P',1)

#Define boundaries
inflow='near(x[0],0)'
outflow='near(x[0],1)'
walls='near(x[1],0)|| near(x[1],1)'

#Define boundary conditions
bcu_noslip = DirichletBC(V,Constant((0,0)),walls)
bcp_inflow = DirichletBC(Q,Constant(8),inflow)
bcp_outflow = DirichletBC(Q,Constant(0),outflow)
bcu = [bcu_noslip]
bcp = [bcp_inflow,bcp_outflow]

#Define trial and test functions
u=TrialFunction(V)
v=TestFunction(V)
p=TrialFunction(Q)
q=TestFunction(Q)

#Define functions for solutions at previous and current time steps
u_n=Function(V)
u_=Function(V)
p_n=Function(Q)
p_=Function(Q)

#Define expressions used in variational forms
U=0.5*(u_n+u)
n=FacetNormal(mesh)
f=Constant((0,0))
k=Constant(dt)
mu=Constant(mu)
rho=Constant(rho)

#Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

#Define stress tensor
def sigma(u,p):
    return 2*mu*epsilon(u)-p*Identity(len(u))

#Define variational problem for step1
F1 = rho*dot((u - u_n) / k, v)*dx + \
     rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
   + inner(sigma(U, p_n), epsilon(v))*dx \
   + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
   - dot(f, v)*dx
a1=lhs(F1)
L1=rhs(F1)

#Define variational problem for step2
a2=dot(nabla_grad(p),nabla_grad(q))*dx
L2=dot(nabla_grad(p_n),nabla_grad(q))*dx-(1/k)*div(u_)*q*dx

#Define variational problem for step3
a3=dot(u,v)*dx
L3=dot(u_,v)*dx-k*dot(nabla_grad(p_-p_n),v)*dx

#Assemble matrices
A1=assemble(a1)
A2=assemble(a2)
A3=assemble(a3)

#Apply boundary conditions to matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

#Create progress bar
#progress=Progress('Time-stepping')
ufile = File("n-s-channel/u.pvd")

#Time-stepping
t=0
for n in range(num_steps):

    #Update current time
    t+=dt
   #Step1:Tentative velocity step
    b1=assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(),b1)

   #Step2:Pressure correction step
    b2=assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(),b2)

   #Step3:Velocity correction step
    b3=assemble(L3)
    solve(A3, u_.vector(),b3)

    #Plot solution
    plot(u_,title='Velocity')
    ufile<<u_

    # Compute error
    u_e = Expression(('4*x[1]*(1.0 - x[1])', '0'), degree=2)
    u_e = interpolate(u_e, V)
    error = (np.abs(u_e.vector() - u_.vector())).max()
    print('t = %.2f: error = %.3g' % (t, error))
    print('max u:', u_.vector().max())

    #Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)