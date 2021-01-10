__author__ = 'Nick_Vasios'
__copyright__ = 'Copyright 2015, All rights reserved'
__email__ = 'vasios@g.harvard.edu'
#
from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# Define tolerance
tol=1.e-12

# Input of user defined parameters
th0=pi/3

# Open File
iout=open('./FDNewton.dat','w')


# Define the dimensions of 'load' and 'displacement' vectors
ndof=1

# Iq is the force distribution vector (needs to be defined explicitly)
iq=np.zeros(ndof)
iq[0]=1.

# a is the dimensionless ``displacement'' vector (no need to define)
a=np.zeros(ndof)

# f is the system of equations (They need to be defined explicitly) [See Function fcn]
f=np.zeros(ndof)

# df is the tangent matrix to the system of equations (Contains derivatives)
# dfinv is the inverse of df
# df elements need to be defined explicitly in function dfcn
df=np.zeros((ndof,ndof))
dfinv=np.zeros((ndof,ndof))

# dls will store the 2 solutions from the 2nd order polynomial w.r.t. ddl
dls=np.zeros(2)

# dao is an araay that stores the last converged ``displacement correction''
dao=np.zeros(ndof)
# al is the dimensionless ``load'' vector
al=0.0

#---------------------------------------------------------------------------------------
# Define the b function needed for calculations
def b(x,y):
	b=1.+x**2.0-2.0*x*sin(y)
	return b
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# Define the system of equations in the form F(u)=0
def fcn(x,y,z):
	bb=b(x[0],y)
	f[0]=(1./sqrt(bb)-1.0)*(sin(y)-x[0])-z
	return f

#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# Define the tangent matrix (stiffness matrix)
# It contains the derivatives of the equations w.r.t the variables
# The function returns both the matrix as df as well as it's inverse as dfinv
def dfcn(x,y,z):
	bb=b(x[0],y)
	# Tangent Matrix
	df[0,0]=1-(1.-sin(y)**2.0)/(bb**1.5)

# Inverse of Tangent Matrix
	dfinv=np.linalg.inv(df)
	return df,dfinv
#---------------------------------------------------------------------------------------


# Define the maximum number of Riks increments
newton=10000
maxiter=100


for i in range(newton):
	if  a[0]>=2.5:
		break
	# Increment starts; Set all variations=0
	da=np.zeros(ndof)
	dl=1.5e-2

	a+=da
	al+=dl
	f=fcn((a),th0,(al))

	fcheck=sqrt(np.dot(f,f))

	if fcheck<tol:
		iout.write(str(a[0])+' '+str(al)+"\n")
		iloop=0
		
	else:
		iters=0
		while fcheck>tol:
			
			iters+=1

			df,dfinv=dfcn((a),th0,(al))
			da=-1.*np.dot(dfinv,f)
			a+=da
			f=fcn((a),th0,(al))
			fcheck=sqrt(np.dot(f,f))
			print a



			if iters>maxiter:
				print 'Convergence cannot achieved within',maxiter,'iterations'
				print 'Program stops'
				exit()
		else:
			iout.write(str(a[0])+' '+str(al)+"\n")
	print 'Newton increment',i,'completed successfully'
print 'The program completed successfully'
iout.close()

with open('./FDNewton.dat') as iread:
    array= np.array([[float(x) for x in line.split()] for line in iread])

exx=np.linspace(0,2.5,501)
exy=np.array([(1./sqrt(b(ax,th0)) - 1.)*(sin(th0)-ax) for ax in exx])


pp=PdfPages('Newton.pdf')
plt.figure()

# Setting text fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=15)

plt.plot(array[:,0],array[:,1],'bo',label='Newton',zorder=1)
plt.plot(exx,exy,label='Exact',linewidth=2,color='black',zorder=0)

# Add some labels
plt.xlabel('$\\alpha$',fontsize=18)
plt.ylabel('$\\lambda$',fontsize=18,rotation=0)

# Add a legend
#axes = plt.gca()
#axes.set_xlim([0,7.5])
#axes.set_ylim([-310,310])

plt.title('$\\theta_0=\\pi/6$') 
plt.legend(loc='best')
pp.savefig()
pp.close()