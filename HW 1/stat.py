import numpy as np 
import matplotlib.pyplot as plt

# A = np.loadtxt("AtomicScale_HW1_data3.txt")
# print len(A)
# plt.plot(A)
# plt.xlabel('UINT')
# #plt.axis([0,900,-9.0,-7.0])
# plt.savefig("Total.png")
# #plt.show()
# A4 = np.loadtxt("AtomicScale_HW1_data4.txt")
# A3 = np.loadtxt("AtomicScale_HW1_data3.txt")
# A2 = np.loadtxt("AtomicScale_HW1_data2.txt")
# A1 = np.loadtxt("AtomicScale_HW1_data1.txt")

# print("4: %d,\n 3: %d,\n 2: %d, \n 1: %d", len(A4), len(A3), len(A2), len(A1))
X = [0.2,0.4,0.6,0.8,1]

Y4 = [1,0.1741,0.30584,0.2648,0.27]
Ys4 = [0.462,0.4569,1,0.9039,0.8675]

Y1 = [0.9981,0.99,1,0.996,0.9947]
Ys1 = [0.9902,0.98376,1,0.99655,0.98696]

Y2 = [1,0.9953,0.9948,0.9934,0.992]
Ys2 = [1,0.983,0.9745,0.9745,0.97355]


plt.plot(X, Y1, 'b-', X, Y2, 'g-')
plt.ylabel("Mean")
plt.title('Convergence of mean')
plt.legend()
plt.savefig("Means1&2.png")
plt.show()

plt.plot(X, Ys1, 'b-', X, Ys2, 'g-')
plt.ylabel("Standard Dev")
plt.title('Convergence of std')
plt.legend()
plt.savefig("STD1&2.png")
plt.show()