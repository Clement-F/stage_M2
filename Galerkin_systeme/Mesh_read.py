import numpy as np
import matplotlib.pyplot as plt

with open('file_data.txt') as f:
    lines = f.readlines()
nb_var = int(lines[0][10:16])
orderx = int(lines[1][10:16])
ordert = int(lines[2][10:16])
nb_cell =int(lines[3][10:16]) ; nb_sub= int(lines[4][13:19])+1
nx = nb_cell * nb_sub
    
sm = int(lines[5][5:11]) -2
with open('mesh_out.txt') as f:
    lines = f.readlines()
    
theta_ = np.zeros((sm,nx))
pos = np.zeros((sm,nx)); LMP = np.zeros((sm,nx)); ent = np.zeros((sm,nx))
ext= np.zeros((sm,nx))
X   = np.zeros(nx)
T   = np.zeros(sm)


for k in range(0,sm):
    print(k,k*(nx+1) )
    print (lines[k*(nx+1)])
    print (lines[k*(nx+1)+1][49])
    T[k] = lines[k*(nx+1)][9:19]
    for i in range(0,(nx)):
        line = lines[k*(nx+1)+i+1]
        X[i] = line[0:10]
        theta_[k,i] = line[10:20]
        pos[k,i] = line[20:30]
        LMP[k,i] = line[29:39]
        ent[k,i] = float(line[39:49])
        if(line[49] == "T") : ext[k,i] =1
        
    
    # plt.plot(X,theta_[k],'b-')
    # plt.show()

xL=-1.; xR=1.
    
plt.pcolormesh(X, T, theta_, shading='auto', cmap='viridis')
plt.xlim(xL,xR)
plt.colorbar()
plt.show()

plt.pcolormesh(X, T, pos, shading='auto', cmap='viridis')
plt.xlim(xL,xR)
plt.colorbar()
plt.show()

plt.pcolormesh(X, T, LMP, shading='auto', cmap='viridis')
plt.xlim(xL,xR)
plt.colorbar()
plt.show()

plt.pcolormesh(X, T, ent, shading='auto', cmap='viridis')
plt.xlim(xL,xR)
plt.colorbar()
plt.show()

plt.pcolormesh(X, T, ext, shading='auto', cmap='viridis')
plt.xlim(xL,xR)
plt.colorbar()
plt.show()
#plt.savefig("mesh_out.png")
#plt.show()