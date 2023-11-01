import numpy as np
import matplotlib.pyplot as plt
h=0.04
k=0.001
x=np.arange(0,2+h,h)
t=np.arange(0,0.045+k,k)
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
# print(T)
matrix1=[]

matrix1.append(T[:,int(0.005/k)])
matrix1.append(T[:,int(0.01/k)])
matrix1.append(T[:,int(0.02/k)])
matrix1.append(T[:,int(0.03/k)])
matrix1.append(T[:,int(0.045/k)])
labels=[0.005,0.01,0.02,0.03,0.045]
print(matrix1)

for i in range(5):
    plt.plot(x,matrix1[i])
    
plt.xlabel('x')
plt.ylabel('phi(x,t)')
plt.legend([f't={value}'for value in labels])
plt.title("explicit method at dt=0.001")
plt.savefig("explicit_method_0.001.png")
plt.show()


    