import numpy as np
import matplotlib.pyplot as plt

with open('file_data.txt') as f:
    lines = f.readlines()

nt = int(lines[0][6:11])
nx = int(lines[1][6:11])
sm = int(lines[2][12:16])

T=[]
for k in range(3,3+sm):
    T.append(float(lines[k][14:22]))

with open('file_sol.txt') as f:
    lines = f.readlines()

U_t = np.zeros((sm,nx))
X   = np.zeros(nx)
dec = 0

def sol_exacte(x,t):
    if(x<0.3-t/2): return 0
    elif(x<0.7-t): return -1
    a = (x-0.7)/t
    if( (-1<a) and (a<0.5)): return a 
    else : return 0.5

U = np.zeros((sm,nx))

for k in range(sm):
    print(k)
    for i in range(1,(nx)):
        X[i] = lines[k*(nx+1) + i][1:10]
        U_t[k][i] = lines[k*(nx+1) +i][10:21]
        
        U[k][i] = sol_exacte(X[i],T[k])
        
        
        # erreur pour burgers 
        
        
    
    plt.plot(X,np.abs(U_t[k]- U[k]),'r')
    plt.plot(X,U[k],'g')
    plt.plot(X,U_t[k],'b')
    plt.show()