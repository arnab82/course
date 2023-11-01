import numpy as np
import matplotlib.pyplot as plt
h=0.04
k=0.001
x=np.arange(0,2,h)
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
A=np.diag([2+2*factor]*(n-2),0)+np.diag([-factor]*(n-3),-1)+np.diag([-factor]*(n-3),+1)
B=np.diag([2-2*factor]*(n-2),0)+np.diag([factor]*(n-3),-1)+np.diag([factor]*(n-3),+1)
for j in range(0,m-1):
    b=T[1:-1,j].copy()
    b=np.dot(B,b)
    b[0]=b[0]+factor*(T[0,j]+T[0,j+1])
    b[-1]=b[-1]+factor*(T[-1,j]+T[-1,j+1])
    solution=np.linalg.solve(A,b)
    # print(solution)
    T[1:-1,j+1]=solution

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
plt.title("Crank Nicholson method at dt=0.001")
plt.savefig("crank_nicholson_method_0.001.png")
plt.show()
x1= np.linspace(0, 2 , 50)
t1 = 0.045
# Calculate Ψ(x, t) using the provided equation
def psi_analytical(x, t):
    result = 0.0
    for n in range(0, len(x)): 
        result += (-1)**n / ((2 * n + 1) ** 2 * (np.pi ** 2)) * np.exp(
                -((2 * n + 1) ** 2 * (np.pi ** 2) * t) / 4) * np.sin((n + 1 / 2) * np.pi * x)
    return result * 8.0
# Calculate Ψ(x, t) for the given x and t
psi_values = psi_analytical(x1, t1)
numerical_values = matrix1[4]

fig, ax1 = plt.subplots()
l1, = ax1.plot(x, psi_values, color='red', label="Analytical")
ax2 = ax1.twinx()
l2, = ax2.plot(x, numerical_values, color='orange', label="Numerical")

# Combine the legends from both axes
lines = [l1, l2]
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels)

ax1.set_xlabel('x')
ax1.set_ylabel('Ψ(x, t)')
ax1.set_title('Wave Function Ψ(x, t)')
ax1.grid(True)

plt.savefig("analytical_solution_vs_numerical_crank_nicholson.png")
plt.show()
import numpy as np
import matplotlib.pyplot as plt
h=0.04
k=0.001
x=np.arange(0,2,h)
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
# #implicit method
A=np.diag([1+2*factor]*(n-2),0)+np.diag([-factor]*(n-3),-1)+np.diag([-factor]*(n-3),+1)
for j in range(1,m):
    b=T[1:-1,j-1].copy()
    b[0]=b[0]+factor*T[0,j]
    b[-1]=b[-1]+factor*T[-1,j]
    solution=np.linalg.solve(A,b)
    print(solution)
    T[1:-1,j]=solution

matrix2=[]

matrix2.append(T[:,int(0.005/k)])
matrix2.append(T[:,int(0.01/k)])
matrix2.append(T[:,int(0.02/k)])
matrix2.append(T[:,int(0.03/k)])
matrix2.append(T[:,int(0.045/k)])
numerical_values1 =  matrix1[4]
numerical_values2 = matrix2[4]

fig, ax1 = plt.subplots()
l1, = ax1.plot(x, numerical_values1, color='red', label="Numerical_Crank_NIcolson")
ax2 = ax1.twinx()
l2, = ax2.plot(x, numerical_values2, color='orange', label="Numerical_implicit")

# Combine the legends from both axes
lines = [l1, l2]
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels)

ax1.set_xlabel('x')
ax1.set_ylabel('Ψ(x, t)')
ax1.set_title('Wave Function Ψ(x, t)')
ax1.grid(True)

plt.savefig("numericall_solution_implicit_vs_numerical_crank_nicholson.png")
plt.show()