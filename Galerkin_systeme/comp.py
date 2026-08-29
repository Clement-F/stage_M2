import numpy as np
import matplotlib.pyplot as plt
p = lambda U:  2.*(U[2] - 0.5*(U[1])**2 /U[0] )

nb_comp = 3
nb_var = []; 
orderx = np.zeros(nb_comp)  ;ordert =np.zeros(nb_comp);
nb_cell = [] ;nb_sub = [];
nx = []     ;sm = [];

lines = []
data_name   = ['file_data.txt','sol_ex/file_blast.txt','blast/conv_data.txt']
sol_name    = ['solSub.txt','sol_ex/blast.dat','blast/conv_sol.txt']
color       = ['b-','r-','g-']


U_t = []
P   = []
X   = []


for k in range(0,nb_comp,1):
    lines= []
    with open(data_name[k])  as f:
        lines.append(f.readlines())
        
    nb_var.append(int(lines[0][0][10:16]))
    orderx[k] = int(lines[0][1][10:16]);   ordert[k] = int(lines[0][2][10:16])
    nb_cell.append(int(lines[0][3][10:16]));   nb_sub.append( int(lines[0][4][13:19]))
    sm.append( int(lines[0][5][5:11]))
    nx.append( int(nb_cell[k] *(nb_sub[k])))
    
    # sm=100
    
    T=[]
    for i in range(6,6+sm[k]):
       T.append(float(lines[0][i][14:30]))

    
    with open(sol_name[k]) as f:
        lines.append(f.readlines())  

    U_t.append([]); P.append([]); X.append([])   


    for i in range(0,sm[k],1):
        U_t[k].append([]); P[k].append([]); X[k].append([])  
        for j in range(0,nx[k]): 
            U_t[k][i].append([])
            X[k][i].append( float(lines[1][i*(nx[k]+1) + j][0:10]))
            U_t[k][i][j].append(float(lines[1][i*(nx[k]+1) +j][10:27]))
            if(nb_var[k] >1) : 
                U_t[k][i][j].append(float(lines[1][i*(nx[k]+1) +j][28:43]))
            if(nb_var[k] >2) :
                U_t[k][i][j].append(float(lines[1][i*(nx[k]+1) +j][44:58]))
                P[k][i].append( p(U_t[k][i][j]))
                
    
for i in range(0,sm[k],1):
    for k in range(0,nb_comp):
        u_ar = np.zeros((nx[k],nb_var[k]));
        pixar= np.zeros((nx[k]))
        
        u_ar[:,:] = np.array(U_t[k][i])
        pixar[:]  = np.array(X[k][i])
        
        plt.plot(pixar,u_ar[:,0] ,color[k])
        # if(nb_var[k] >1) :
        #     plt.plot(X,U_t[i,:,1] ,'b-')
        # if(nb_var[k] >2) :
        #     plt.plot(X,U_t[i,:,2] ,'b-')
    plt.show()
              
    



if(False):
    with open('sol_ex/buckley.dat') as f:
        lines = f.readlines()
        
    X_B   = np.zeros(30000)
    U_B   = np.zeros(30000)
    for i in range(0,30000):
        X_B[i] = lines[i][0:16]
        U_B[i] = lines[i][17:32]
        
    for i in range(0,nb_cell+1): 
        1+1
        #plt.plot([X_cell[i],X_cell[i]],[-1,2],linestyle='--', color='gray')
    plt.plot(X_B,U_B, 'r')
    
    plt.plot(X,U_t[k],'b', marker=".")
    #plt.plot(X,U_t2[k],'go')
    plt.ylim(-.2,1.2); plt.xlim(-1.,1.)
    #plt.show()     
    plt.savefig("endfig_save.png")    


if(False):
    with open('sol_ex/sod_modif_den.dat') as f:
        lines = f.readlines()
    
    rho_ex = np.zeros(5000)
    
    X_ex   = np.zeros(5000)
    rho_ex = np.zeros(5000)
    Xvel_ex= np.zeros(9998)
    vel_ex = np.zeros(9998)
    pre_ex = np.zeros(5000)
    
    for i in range(0,5000):
        X_ex[i] = lines[i][0:18]
        rho_ex[i] = lines[i][17:32]
    
    with open('sol_ex/sod_modif_vel.dat') as f:
        lines = f.readlines()
        
    for i in range(0,9998):
        Xvel_ex[i] = lines[i][0:17]
        vel_ex[i] = lines[i][17:32]
    
    with open('sol_ex/sod_modif_pre.dat') as f:
        lines = f.readlines()
        
    for i in range(0,5000):
        pre_ex[i] = lines[i][17:32]
        
    plt.plot(X,U_t[k,:,0],'b', marker='.')
    plt.plot(X_ex,rho_ex,'b')
    # plt.plot(Xvel_ex,vel_ex,'r')
    # plt.plot(X,U_t[k,:,1]/U_t[k,:,0],'r', marker='.')
    # plt.plot(X_ex,pre_ex,'k')
    # plt.plot(X,P_[k,:],'k', marker='.')
    
    plt.savefig("endfig_save.png")    

    