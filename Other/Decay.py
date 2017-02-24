import numpy as np 
import sys
import matplotlib.pylab as plt

var = sys.argv
lbda = float(var[1])
if lbda >=1:
	print "Lambda should be smaller than 1! Try again!"
	sys.exit()

Nnuc = int(var[2])
time = 0
N =[]
rate=[]
# print lbda, Nnuc
# print type(lbda), type(Nnuc)
while(Nnuc>0):
	time+=1
	num= Nnuc
	N.append(Nnuc)
	for i in range(num):
		r = np.random.rand()
		if r<lbda:
			Nnuc-=1
	rate.append(num-Nnuc)

t = np.arange(time)

lN = np.log(N) 

plt.plot(t,N)
plt.xlabel("time")
plt.ylabel("No. of particles remaining")
plt.savefig("Nnuc_vs_time.jpg")
plt.show()


plt.plot(t,lN)
plt.xlabel("time")
plt.ylabel("Log of No. of particles remaining")
plt.savefig("Log_Nnuc_vs_time.jpg")
plt.show()


plt.plot(t,rate)
plt.xlabel("time")
plt.ylabel("Rate of Decay")
plt.savefig("Decay_rate_vs_time.jpg")
plt.show()
