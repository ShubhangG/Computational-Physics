import numpy as np 
import matplotlib.pylab as plt

theta = np.random.rand(90000)*2*np.pi


x=0
y=0
x1=y1=0
N=np.sqrt(np.arange(90000.0))
R1=[]
R= []
X=[]
Y=[]

for step in theta:
	del_x1 = np.cos(step)
	del_y1 = np.sin(step)
	del_x = np.random.uniform(-np.sqrt(2),np.sqrt(2))
	del_y = np.random.uniform(-np.sqrt(2),np.sqrt(2))
	x = x + del_x
	y = y + del_y
	x1 = x1 + del_x1
	y1 = y1 + del_y1
	X.append(x)
	Y.append(y)
	R1.append(np.sqrt(x1*x1 + y1*y1))
	R.append(np.sqrt(x*x + y*y))

plt.plot(N,R,'b-',N,R1,'g-',N,N,'r--')
plt.axis([0,300,0,300])
plt.xlabel('Sqrt N')
plt.ylabel('R')
plt.savefig('randomW_RvssqrtN.jpg')
plt.show()
