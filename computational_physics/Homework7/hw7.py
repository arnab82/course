import numpy as np
import matplotlib.pyplot as plt
h=0.045
k=0.0015
x=np.arange(0,2+h,h)
t=np.arange(0,0.05+k,k)
boundary_conditions=[0,0]
n=len(x)
m=len(t)
T=np.zeros((n,m))
T[0,:]=boundary_conditions[0]
T[-1,:]=boundary_conditions[1]
#initial conditions
for i in range(n):
    if 0 <= x[i] <= 1:
        T[i] = x[i]
    else:
        T[i] = -x[i] + 2
print(T)
factor=k/h**2
#explicit method
for j in range(1,m):
    for i in range(1,n-1):
        # print(i,j)
        T[i,j]=factor*T[i-1,j-1]+(1-2*factor)*T[i,j-1]+factor*T[i+1,j-1]
print(T)
R=np.linspace(1,0,m)
B=np.linspace(0,1,m)
G=0
for j in range(m):
    plt.plot(x,T[:,j],color=[R[j],G,B[j]])
plt.xlabel('distance[m]')
plt.ylabel('Temperature[$\degree$ C]')
plt.legend([f't={value}s'for value in t.round(3)])
plt.show()

    