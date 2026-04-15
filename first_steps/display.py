import numpy as np
import matplotlib.pyplot as plt

with open('file_sol.txt') as f:
    lines = f.readlines()

nt = len(lines)
nx = len(lines[0])
U_t = np.zeros((nt,int(nx/12)))

for k in range(0,nt):
    U_t[k,0]= (float(lines[k][1:9]))

    for i in range(9,nx,8+4):
        # print(i,(lines[k][i:i+8+4]) )
        # print(i,float(lines[k][i:i+8+4]) )

        U_t[k,int((i-1)/(8+4))] =(float(lines[k][i:i+8+4]))

X = np.linspace(0,1,int(nx/12))

for i in range(0,nt,10):
    plt.plot(X,U_t[i])
    plt.show()