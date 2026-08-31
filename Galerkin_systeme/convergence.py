import numpy as np
import matplotlib.pyplot as plt


with open('file_data.txt') as f:
    lines = f.readlines()
nb_var = int(lines[0][10:16])
orderx = int(lines[1][10:16])
ordert = int(lines[2][10:16])
nb_cell =int(lines[3][10:16]) ; nb_sub= int(lines[4][13:19])
nx = nb_cell * nb_sub

with open('convergence_err.txt') as f:
    lines = f.readlines()

size = len(lines)
nx =[]; err_L1 =[]; err_L2 =[]; err_Li=[] 
ord1= []; ord2= []; ord3= []; ord4 =[]; ord5 =[]; ord6= []

vit1=0. ; vit2=0.; viti=0.

print('h | err1     | O1    | err2  | O2 \n')

for k in range(0,int(size/7)):
    nx.append(      int(lines[7*k+2][9:15])  )
    err_L1.append(  float(lines[7*k+3][11:30])) 
    err_L2.append(  float(lines[7*k+4][11:30])) 
    err_Li.append(  float(lines[7*k+5][11:30])) 
    ord1.append(  pow(nx[-1],-1))
    ord2.append(  pow(nx[-1],-2))
    ord3.append(  pow(nx[-1],-3))
    ord4.append(  pow(nx[-1],-4))
    if k>0 :
        vit1 = np.log(err_L1[k]/err_L1[k-1])/np.log(nx[k]/nx[k-1])
        vit2 = np.log(err_L2[k]/err_L2[k-1])/np.log(nx[k]/nx[k-1])
        viti = np.log(err_Li[k]/err_Li[k-1])/np.log(nx[k]/nx[k-1])
        
        print(1./nx[k],"|",err_L1[k],"|",vit1,"|",err_L2[k],"|",vit2)
      



plt.loglog(nx,ord1,linestyle ='-.', color='lightgray')
plt.loglog(nx,ord2,linestyle ='-.', color='lightgray')
plt.loglog(nx,ord3,linestyle ='-.', color='lightgray')
plt.loglog(nx,ord4,linestyle ='-.', color='lightgray')
# plt.loglog(nx,err_L1,linestyle ='-',marker='x', label="err_L1")
plt.loglog(nx[0:5],err_L2[0:5],linestyle ='-',marker='x', label="err_L2")
plt.loglog(nx[5:10],err_L2[5:10],linestyle ='-',marker='x', label="LMP")
plt.loglog(nx[10:15],err_L2[10:15],linestyle ='-',marker='x', label="LMP + Relax")
# plt.loglog(nx[18:],err_L2[18:],linestyle ='-',marker='x', label="GMP")
# plt.loglog(nx,err_Li,linestyle ='-',marker='x', label="err_Linf")
plt.ylabel("erreur")
plt.xlabel("DoF")
title = "convergence de schéma monolithique d'ordre "+str(orderx)

plt.title(title)

plt.legend()
plt.grid()
plt.savefig("convergence_plot.png")
#plt.show()
