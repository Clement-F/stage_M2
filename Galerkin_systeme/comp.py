import numpy as np
import matplotlib.pyplot as plt
p = lambda U:  2.*(U[2] - 0.5*(U[1])**2 /U[0] )

nb_comp = 3
nb_var = []; 
orderx = [] ;ordert = []
nb_cell = [] ;nb_sub = [];
nx = []     ;sm = [];

lines = []

data_name   = ['pos_data.txt','LMP_data.txt','LMP_Rel_pos_data.txt']
sol_name    = ['sol_pos.txt','sol_LMP.txt','sol_LMP_Rel_pos.txt']


# data_name   = ['file_data.txt','sol_ex/file_blast.txt','blast/conv_data.txt']
# sol_name    = ['solSub.txt','sol_ex/blast.dat','blast/conv_sol.txt']
color       = ['b','r','g','y','k']
markers     = ['.','.','.','.','None']
linestyle   = ['None','-','-','-','-.']

# legend      = ['ordre 1','ordre 3', 'ordre 9','sol exacte']
legend      = ['pos','PML','PML+Rel+pos','exacte']

U_t = []
P   = []
X   = []


for k in range(0,nb_comp,1):
    print(sol_name[k])
    lines= []
    with open(data_name[k])  as f:
        lines.append(f.readlines())
        
    nb_var.append(int(lines[0][0][10:16]))
    orderx.append(int(lines[0][1][10:16]));    ordert.append( int(lines[0][2][10:16]))
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
              
        plt.plot(pixar,u_ar[:,0] ,color[k], marker=markers[k], markersize=3, linestyle=linestyle[k], label=legend[k])
        # if(nb_var[k] >1) :
        #     plt.plot(X,U_t[i,:,1] ,'b-')
        # if(nb_var[k] >2) :
        #     plt.plot(X,U_t[i,:,2] ,'b-')
        
        # if(k==0): plt.plot(pixar,np.zeros((nx[0])),color='k',linestyle='--')  
    plt.ylim(-.2,1.2)
    plt.legend()
    # plt.axis('off')
    plt.show()
              
    

