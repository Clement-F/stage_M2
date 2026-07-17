import numpy as np
import matplotlib.pyplot as plt

with open('file_data.txt') as f:
    lines = f.readlines()
nb_var = int(lines[0][10:16])
orderx = int(lines[1][10:16])
ordert = int(lines[2][10:16])
nb_cell =int(lines[3][10:16]) ; nb_sub= int(lines[4][13:19])
nx = nb_cell * nb_sub
    
sm = int(lines[5][5:11]) -1

with open('mesh_out.txt') as f:
    lines = f.readlines()
    
theta_ = np.zeros((sm,nx))
X   = np.zeros(nx)
T   = np.zeros(sm)


for k in range(0,sm):
    print(k)
    print (lines[k*(nx+1)])
    T[k] = lines[k*(nx+1)][9:19]
    for i in range(0,(nx)):
        X[i] = lines[k*(nx+1)+i+1][0:10]
        theta_[k,i] = lines[k*(nx+1) +i+1][11:23]
        
    
    # plt.plot(X,theta_[k],'b-')
    # plt.show()
    
    
plt.figure(figsize=(8, 6))
plt.pcolormesh(X, T, theta_, shading='auto', cmap='viridis')
plt.savefig("mesh_out.png")
#plt.show()