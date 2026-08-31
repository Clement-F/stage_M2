import numpy as np
import matplotlib.pyplot as plt
p = lambda U:  2.*(U[2] - 0.5*(U[1])**2 /U[0] )

nb_comp = 1
nb_var = []; 
orderx = [] ;ordert = []
nb_cell = [] ;nb_sub = [];
nx = []     ;sm = [];

lines = []


# data_name   = ['1D_data.txt','3D_data.txt','6D_data.txt','1D_data.txt']
# sol_name    = ['sol1D.txt','sol3D.txt','sol6D.txt','solex.txt']

# data_name   = ['pos_data.txt','LMP_pos_data.txt','LMP_pos_Rel_data.txt','solex_data.txt']
# sol_name    = ['sol_pos.txt','sol_LMP_pos.txt','sol_LMP_pos_Rel.txt','solex.txt']

data_name   = ['blast_data.txt','unif_data.txt','blast_data.txt']
sol_name    = ['blast_file.txt','sol_unif.txt','blast_file.txt']



color       = ['b','r','g','c','k']
markers     = ['.','.','.','x','None']
linestyle   = ['-','-','-','-','-']

# legend      = ['ordre 1','ordre 3', 'ordre 6','sol exacte']
legend      = ['LMP+pos','PML+Rel+pos','PML Rel','sol exacte']

U_t = []
P   = []
X   = []
T   = []
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
    
    # if(k==3): nx[k] =nb_cell[k]*(orderx[k]+1)
    
    # sm=100
    if(k==0):
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
    print (i)
    for k in range(0,nb_comp):
        u_ar = np.zeros((nx[k],nb_var[k]));
        p_ar = np.zeros((nx[k]))
        pixar= np.zeros((nx[k]))
        
        u_ar[:,:] = np.array(U_t[k][i])
        p_ar[:] = np.array(P[k][i])
        pixar[:]  = np.array(X[k][i])
              
        plt.plot(pixar,u_ar[:,0] ,color[k], marker=markers[k], markersize=3, linestyle=linestyle[k], label=legend[k])
        # plt.plot(pixar,p_ar[:] ,color[k], marker=markers[k], markersize=3, linestyle=linestyle[k], label=legend[k])
        # if(nb_var[k] >1) :
        #     plt.plot(X,U_t[i,:,1] ,'b-')
        # if(nb_var[k] >2) :
        #     plt.plot(X,U_t[i,:,2] ,'b-')
        
        # if(k==0): plt.plot(pixar,np.zeros((nx[0])),color='k',linestyle='--')  
    # plt.ylim(3.8,4.8)
    # plt.xlim(0.2,0.8)
    plt.legend(loc=0)
    plt.title("t = "+str(round(T[i],3) ))
    # plt.axis('off')
    plt.show()
              
    

