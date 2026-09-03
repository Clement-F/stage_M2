import numpy as np
import matplotlib.pyplot as plt

with open('file_data.txt') as f:
    lines = f.readlines()
nb_var = int(lines[0][10:16])
orderx = int(lines[1][10:16])
ordert = int(lines[2][10:16])
nb_cell =int(lines[3][10:16]) ; nb_sub= int(lines[4][13:19])
nx = nb_cell * (nb_sub+1)
    
#sm = 40
sm = int(lines[5][5:11]) -2
with open('mesh_out.txt') as f:
    lines = f.readlines()
    
theta_ = np.zeros((sm,nx,nb_var))
pos = np.zeros((sm,nx,nb_var)); LMP = np.zeros((sm,nx,nb_var))
ext= np.zeros((sm,nx))
X   = np.zeros(nx)
T   = np.zeros(sm)
line =[]
p=0

for k in range(0,sm):
    T[k] = lines[k*(nx+1)][9:19]
    for i in range(0,(nx)):
        # print(p)
        line = lines[k*(nx+1) +i+1]
        X[i] = lines[k*(nx+1)+i+1][0:8]
        p=8 
        for j in range(nb_var):
            # print(p)
            theta_[k,i,j] = float(line[p:p+6])
            p +=6
        p+=1
        for j in range(nb_var):
            # print(p)
            if(line[p:p+6] != '***** ') :pos[k,i,j] = float(line[p:p+6])
            p +=6
        p+=1
            
        for j in range(nb_var):
            # print(p)
            LMP[k,i,j] = float(line[p:p+6])
            p +=6

        
        if(lines[k*(nx+1)+i+1][p+1] == "T") : ext[k,i] =1
        
    
    # plt.plot(X,theta_[k],'b-')
    # plt.show()

xL=-1.; xR=1.

for i in range(0,nb_var):
    plt.pcolormesh(X, T, theta_[:,:,i], shading='auto', cmap='viridis')
    plt.xlim(xL,xR)
    plt.title("theta")
    plt.colorbar()
    plt.show()

    

    plt.pcolormesh(X, T, pos[:,:,i], shading='auto', cmap='viridis')
    plt.xlim(xL,xR)
    plt.title("theta pos")
    plt.colorbar()
    plt.show()
    
    plt.pcolormesh(X, T, LMP[:,:,i], shading='auto', cmap='viridis')
    plt.xlim(xL,xR)
    plt.title("theta PML")
    plt.colorbar()
    plt.show()
    
    plt.pcolormesh(X, T, ext, shading='auto', cmap='viridis')
    plt.xlim(xL,xR)
    plt.title("theta relax")
    plt.colorbar()
    plt.show()
#plt.savefig("mesh_out.png")
#plt.show()
