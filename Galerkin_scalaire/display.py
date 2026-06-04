import numpy as np
import matplotlib.pyplot as plt

with open('file_data.txt') as f:
    lines = f.readlines()
order = int(lines[0][10:16])
nx = int(lines[2][5:11])
sm = int(lines[3][5:11])

T=[]
#for k in range(4,4+sm):
 #   T.append(float(lines[k][14:30]))


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
        X[i] = lines[k*(nx+1) + i][0:10]
        U_t[k][i]   = lines[k*(nx+1) +i][11:23]
        U_ex[k][i]  = lines[k*(nx+1) +i][27:43]        
        err[k][i] = abs(U_t[k][i] - U_ex[k][i]) 
                   
        
    err_L2[k] = np.sqrt(sum(err[k][:]**2)) * 1/nx
    #plt.plot(X,err[k],'r')
    plt.plot(X,U_ex[k],'g')
    plt.plot(X,U_t[k],'b')
    plt.ylim(-1.2,1.2)
    plt.show()         
        
#for k in range(1,sm):
    #err_L2[k] = np.sqrt(sum(err[k][:]**2)) * 1/nx
    #plt.plot(X,err[1],'r')
#plt.plot(X,U_ex[0],'g')
#plt.plot(X,U_t[0],'k')
#plt.plot(X,U_t[-1],'b')
#plt.plot(X,U_t[2],'r')
#plt.xlim(0.0,0.1)
#print(sum(U_t[0]), sum(U_t[-1]))
#print(max(err_L2))
    