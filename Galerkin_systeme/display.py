import numpy as np
import matplotlib.pyplot as plt


p = lambda U:  0.4*(U[2] - 0.5*(U[1])**2 /U[0] )
    
with open('file_data.txt') as f:
    lines = f.readlines()
    
nb_var = int(lines[0][10:16])
orderx = int(lines[1][10:16])
ordert = int(lines[2][10:16])
nb_cell =int(lines[3][10:16]) ; nb_sub= int(lines[4][13:19])
nx = nb_cell *nb_sub

sm = int(lines[5][5:11])

T=[]
for k in range(6,6+sm):
   T.append(float(lines[k][14:30]))


#   solDG   SolSub   file_sol    

with open('solSub.txt') as f:
    lines1 = f.readlines()   
with open('solSub'+'ex'+'.txt') as f:
    lines2 = f.readlines()   


U_t = np.zeros((sm,nx,nb_var));     P_ = np.zeros((sm,nx))
U_ex = np.zeros((sm,nx,nb_var));    P2_ = np.zeros((sm,nx))
X   = np.zeros(nx)

X_cell = np.zeros(nb_cell+1);       U_cell=  np.zeros((sm,nb_cell,nb_var))
X_midcell = np.zeros(nb_cell)
dec = 0

X_cell[-1]=1.


for k in range(0,sm,1):
    for i in range(0,nx): 
        X[i] = lines1[k*(nx+1) + i][0:10]
        U_t[k][i][0]   = lines1[k*(nx+1) +i][10:27]
        if(nb_var >1) : 
            U_t[k][i][1]   = lines1[k*(nx+1) +i][28:43]
        if(nb_var >2) :
            U_t[k][i][2]   = lines1[k*(nx+1) +i][44:58]
            
        # if(i%nb_sub==0) : 
        #     X_cell[int(i/nb_sub)] =lines1[k*(nx+1) + i][59:69]
        #     U_cell[k][int(i/nb_sub)-1][0] = sum(U_t[k,i-nb_sub:i,0])/nb_sub
        #     X_midcell[int(i/nb_sub)] = (X_cell[int(i/nb_sub)]+X_cell[int(i/nb_sub)+1])/2
            
    
        #U_ex[k][i][0]   = lines2[k*(nx+1) +i][11:27]
        # if(nb_var >1) : 
        #      U_ex[k][i][1]   = lines2[k*(nx+1) +i][28:43]
        # if(nb_var >2) :
        #     U_ex[k][i][2]   = lines2[k*(nx+1) +i][44:58]
            
        
        #P_[k,i] = p(U_t[k,i,:])
        #U_cell[k][-1][0] = sum(U_t[k,-nb_sub-1:-1,0])/nb_sub
    #U_cell[k][-1][0] = sum(U_t[k,-nb_sub-1:-1,0])/nb_sub
                  
    for i in range(0,nb_cell+1): 
        1+1
        #plt.plot([X_cell[i],X_cell[i]],[-1,2],linestyle='--', color='gray')
    #plt.plot(X_midcell,U_cell[k,:,0], marker='.')
    #plt.plot(X,U_t[k,:,0],'b-.')
    plt.plot(X,U_t[k,:,0],'b', marker='.')
    plt.plot(X,U_ex[k,:,0],'g-')
    #plt.plot(X,abs(U_t[k,:,0]-U_ex[k,:,0]),'k-')
    #print(max(abs(U_t[k,:,0]-U_ex[k,:,0])))
    #print(sum(U_t[k,:,0]))
    if(nb_var>1):
        1+1
        #plt.plot(X,U_ex[k,:,1],'k-')
        #plt.plot(X,U_t[k,:,1],'r-', marker='.')
        #plt.plot(X,abs(U_t[k,:,1]-U_ex[k,:,1]),'k-')
        #plt.plot(X,U_t[k,:,1]/U_t[k,:,0],'r-', marker='.')
    if(nb_var>2):
        1+1
        #plt.plot(X,U_ex[k,:,2],'k-')
        #plt.plot(X,U_t[k,:,2],'g-', marker='.')
        #plt.plot(X,abs(U_t[k,:,2]-U_ex[k,:,2]),'k-')
    
    #plt.plot(X,P_[k,:],'k-x')
    #m=-.2; M=1.2
    
    #plt.ylim(3.8,4.8); plt.xlim(0.2,0.8)
   # plt.ylim(m,M);plt.xlim(-0.,1.)
    #plt.ylim(-.2, 7.2)
    plt.show()    
    # plt.savefig("save_"+str(T[k])+".png")    
    # plt.cla() 
    
# plt.plot(X,U_t[k,:,0],'b', marker='.')
# plt.ylim(3.,4.8);plt.xlim(0.5,2.3)
# plt.savefig("save_"+str(T[k])+"zoom.png")    
    



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
    plt.plot(Xvel_ex,vel_ex,'r')
    plt.plot(X,U_t[k,:,1]/U_t[k,:,0],'r', marker='.')
    plt.plot(X_ex,pre_ex,'k')
    plt.plot(X,P_[k,:],'k', marker='.')
    
    plt.savefig("endfig_save.png")    

    