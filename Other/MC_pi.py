import numpy as np 
import sys
import matplotlib.pylab as plt

N = 900
N_sp = 0.0
r = np.random.rand(N)
A = []
p=[]
for i in range(N):
	x_i= np.random.uniform(-1,1)
	y_i= np.random.uniform(-1,1)
	p.append(np.pi)
	if x_i*x_i + y_i*y_i < 1:
		N_sp+=1
	if i==0:
		A.append(0)
		continue
	A.append(4*N_sp/i)
	
print A[i]
n = np.arange(N)
plt.plot(n,A,'b-',n,p,'r--')
plt.xlabel("no. of iterations")
plt.ylabel("Area value")
plt.savefig("Pi_determination.png")
plt.show()

#Area = N_sp/N*4

#print Area