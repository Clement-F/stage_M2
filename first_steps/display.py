import numpy as np
import matplotlib.pyplot as plt

with open('file_data.txt') as f:
    lines = f.readlines()

nt = int(lines[0][5:11])
nx = int(lines[1][5:11])
sm = int(lines[2][12:16])

T=[]
for k in range(3,3+sm):
    T.append(float(lines[k][14:22]))

with open('file_sol.txt') as f:
    lines = f.readlines()

U_t = np.zeros((sm,nx))
X   = np.zeros(nx)
dec = 0

U_ex = np.zeros((sm,nx))
err = np.zeros((sm,nx))
err_L2 = np.zeros(sm)

for k in range(0,sm):
    for i in range(0,(nx)):
        X[i] = lines[k*(nx+1) + i][1:11]
        U_t[k][i]   = lines[k*(nx+1) +i][11:23]
        U_ex[k][i]     = lines[k*(nx+1) +i][23:35] 
        
        err[k][i] = abs(U_t[k][i] - U_ex[k][i])
                        
        
    err_L2[k] = np.sqrt(sum(err[k][:]**2)) * 1/nx
    plt.plot(X,err[k],'r')
    #plt.plot(X,U_ex[k],'g')
    #plt.plot(X,U_t[k],'b')
    #plt.xlim(0.0,0.1)
    plt.show()

print(max(err_L2))
    