__author__ = 'Nick_Vasios'
__copyright__ = 'Copyright 2015, All rights reserved'
__email__ = 'vasios@g.harvard.edu'
#
from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# Define tolerance
tol=1.e-9

# Input of user defined parameters
th0=pi/3
w=0.25
psi=1.0
dll=2.5e-4

# Open File
iout=open('./FD_'+str(w)+'.dat','w')


# Define the dimensions of 'load' and 'displacement' vectors
ndof=2

# Iq is the force distribution vector (needs to be defined explicitly)
iq=np.zeros(ndof)
iq[0]=0.
iq[1]=1.

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
def fcn(x,y,z,w):
	bb=b(x[0],y)
	f[0]=-w*(x[1]-x[0])+(1./sqrt(bb)-1.0)*(sin(y)-x[0])
	f[1]=w*(x[1]-x[0])-z
	return f

#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# Define the tangent matrix (stiffness matrix)
# It contains the derivatives of the equations w.r.t the variables
# The function returns both the matrix as df as well as it's inverse as dfinv
def dfcn(x,y,z,w):
	bb=b(x[0],y)
	# Tangent Matrix
	df[0,0]=(1.+w)-(1.-sin(y)**2.0)/(bb**1.5)
	df[0,1]=-w
	df[1,0]=-w
	df[1,1]=w
# Inverse of Tangent Matrix
	dfinv=np.linalg.inv(df)
	return df,dfinv
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# The arc length function solves the 2nd order equation w.r.t. ddl
# Returns two values as ddl1 and ddl2    
def arc(da,dab,dat,dl,iq):

	#Arc Length Parameters
	psi=1.0
	dll=1.e-3

	# Calculate the coefficients of the polynomial
	c1=np.dot(dat,dat)+psi**2.0*np.dot(iq,iq)
	c2=2.*(np.dot((da+dab),dat)+dl*psi**2.0*np.dot(iq,iq))
	c3=np.dot((da+dab),(da+dab))+dl**2.0*psi**2.0*np.dot(iq,iq)-dll**2.0

	if c2**2.0-4.*c1*c3>0.:
		dls=np.roots([c1,c2,c3])
		ddl1=dls[0]
		ddl2=dls[1]
	else:
		ddl1=-c2/2*c1
		ddl2=-c2/2*c1
		print 'Possible issue in Arc Length equation'

	return ddl1,ddl2
#---------------------------------------------------------------------------------------


# Define the maximum number of Riks increments
riks=20000
maxiter=10


for i in range(riks):
	if  a[1]>=3.5:
		break
	# Increment starts; Set all variations=0
	da=np.zeros(ndof)
	dab=np.zeros(ndof)
	dat=np.zeros(ndof)
	dda1=np.zeros(ndof)
	dda2=np.zeros(ndof)
	dda=np.zeros(ndof)
	dalfa=np.zeros(ndof)
	dl=0.


	df,dfinv=dfcn((a+da),th0,(al+dl),w)
	dat=np.dot(dfinv,iq)
	
	ddl1,ddl2=arc(da,dab,dat,dl,iq)

	dda1=dab+ddl1*dat
	dda2=dab+ddl2*dat
	
	det=np.linalg.det(df)
	
	if np.sign(det)==np.sign(ddl1):
		dda=dda1
		ddl=ddl1
	else:
		dda=dda2
		ddl=ddl2

	dalfa=da+dda
	dlamda=dl+ddl

	f=fcn((a+dalfa),th0,(al+dlamda),w)

	fcheck=sqrt(np.dot(f,f))

	if fcheck<tol:
		a=a+dalfa
		al=al+dlamda
		iout.write(str(a[1])+' '+str(al)+"\n")
		dao=dalfa
		dlo=dlamda
		iloop=0
		
	else:
		iters=0
		while fcheck>tol:
			
			iters+=1
			da=dalfa
			dl=dlamda

			f=fcn((a+da),th0,(al+dl),w)
			df,dfinv=dfcn((a+da),th0,(al+dl),w)

			dab=-np.dot(dfinv,f)
			dat=np.dot(dfinv,iq)

			ddl1,ddl2=arc(da,dab,dat,dl,iq)

			dda1=dab+ddl1*dat
			dda2=dab+ddl2*dat

			det=np.linalg.det(df)

			daomag=np.dot(dao,dao)

			if daomag==0.:
				if np.sign(dl+ddl1)==np.sign(det):
					dda=dda1
					ddl=ddl1
				else:
					dda=dda2
					ddl=ddl2
			else:
				aux1=np.dot(da+dda1,dao)
				aux2=np.dot(da+dda2,dao)

				aux3=dlamda*(dl+ddl1)*np.dot(iq,iq)
				aux4=dlamda*(dl+ddl2)*np.dot(iq,iq)

				dot1=aux1+psi**2.*aux3
				dot2=aux2+psi**2.*aux4

				if dot1>dot2:
					dda=dda1
					ddl=ddl1
				else:
					dda=dda2
					ddl=ddl2

			if ddl1==ddl2:
				dda=dda1
				ddl=ddl1

			dalfa=da+dda
			dlamda=dl+ddl

			f=fcn((a+dalfa),th0,(al+dlamda),w)

			fcheck=np.linalg.norm(f)

			if iters>maxiter:
				iters=maxiter+1
				break

		if iters>maxiter:
			print 'Convergence cannot achieved within',maxiter,'iterations'
			print 'Program stops'
			exit()
		else:
			a+=dalfa
			al+=dlamda
			iout.write(str(a[1])+' '+str(al)+"\n")
			dao=dalfa
			dlo=dlamda
	#print 'Riks increment',i,'completed successfully'
print 'The program completed successfully'
iout.close()
 
#iread=open('Riks.dat', 'r')
with open('./FD_'+str(w)+'.dat') as iread:
    u, l = [float(x) for x in iread.readline().split()]
    array= np.array([[float(x) for x in line.split()] for line in iread])


plt.figure()

# Setting text fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=15)

plt.plot(array[:,0],array[:,1],label='$\\beta=0.3$',linewidth=2)

# Add some labels
plt.xlabel('$\\alpha$',fontsize=18)
plt.ylabel('$\\lambda$',fontsize=18,rotation=0)
plt.title('$\\theta_0=\\pi/6$') 
plt.legend(loc='best')
plt.show()