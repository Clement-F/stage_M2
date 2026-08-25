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
    
theta_ = np.zeros((sm,nx,3))
pos = np.zeros((sm,nx,3)); LMP = np.zeros((sm,nx,3))
ext= np.zeros((sm,nx))
X   = np.zeros(nx)
T   = np.zeros(sm)
line =[]

for k in range(0,sm):
    print(k,k*(nx+1) )
    # print (lines[k*(nx+1)])
    #print (lines[k*(nx+1)+0+1][65])
    T[k] = lines[k*(nx+1)][9:19]
    for i in range(0,(nx)):
        line = lines[k*(nx+1) +i+1]
        X[i] = lines[k*(nx+1)+i+1][0:8]
        theta_[k,i,0] = float(line[8:14])
        theta_[k,i,1] = float(line[14:20])
        theta_[k,i,2] = float(line[20:26])
        
        if(line[27:33] != "***** ") :pos[k,i,0] = float(line[27:33])
        if(line[33:39] != "***** ") :pos[k,i,1] = float(line[33:39])
        if(line[39:45] != "***** ") :pos[k,i,2] = float(line[39:45])
        
        LMP[k,i,0] = float(line[46:52])
        LMP[k,i,1] = float(line[52:58])
        LMP[k,i,2] = float(line[58:64])
        
        if(lines[k*(nx+1)+i+1][65] == "T") : ext[k,i] =1
        
    
    # plt.plot(X,theta_[k],'b-')
    # plt.show()

xL=-0.; xR=1.

for i in range(0,3):
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
    plt.title("theta LMP")
    plt.colorbar()
    plt.show()
    
    plt.pcolormesh(X, T, ext, shading='auto', cmap='viridis')
    plt.xlim(xL,xR)
    plt.title("theta relax")
    plt.colorbar()
    plt.show()
#plt.savefig("mesh_out.png")
#plt.show()