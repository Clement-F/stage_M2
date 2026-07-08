import numpy as np
import matplotlib.pyplot as plt


p = lambda U:  2.*(U[2] - 0.5*(U[1])**2 /U[0] )
    
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

with open('solSub'+'ex'+'.txt') as f:
    lines2 = f.readlines()   

with open('solSub.txt') as f:
    lines1 = f.readlines()   

U_t = np.zeros((sm,nx,nb_var)); P_ = np.zeros((sm,nx))
U_2t = np.zeros((sm,nx,nb_var)); P2_ = np.zeros((sm,nx))
X   = np.zeros(nx)
dec = 0




for k in range(0,sm):
    for i in range(0,nx): 
        X[i] = lines1[k*(nx+1) + i][0:10]
        U_t[k][i][0]   = lines1[k*(nx+1) +i][11:27]
        if(nb_var >1) : 
            U_t[k][i][1]   = lines1[k*(nx+1) +i][28:43]
        if(nb_var >2) :
            U_t[k][i][2]   = lines1[k*(nx+1) +i][44:58]
            
            
        U_2t[k][i][0]   = lines2[k*(nx+1) +i][11:27]
        # print(max(abs(U_t[k,:,0]-U_2t[k,:,0])))
        if(nb_var >1) : 
             U_2t[k][i][1]   = lines2[k*(nx+1) +i][28:43]
        #     print(max(abs(U_t[k,:,1]-U_2t[k,:,1])))
        if(nb_var >2) :
            U_2t[k][i][2]   = lines2[k*(nx+1) +i][44:58]
        #     print(max(abs(U_t[k,:,2]-U_2t[k,:,2])))
            
        
        #P_[k,i] = p(U_t[k,i,:])
                  
        
    plt.plot(X,U_t[k,:,0],'b', marker='.')
    #plt.plot(X,-np.sin(X),'g')
    #plt.plot(X,U_t[0,:,0],'g-')
    plt.plot(X,U_2t[k,:,0],'g-')
    #plt.plot(X,abs(U_t[k,:,0]-U_2t[k,:,0]),'k-')
    if(nb_var>1):
        1+1
        plt.plot(X,U_t[k,:,1],'r-')
        plt.plot(X,U_2t[k,:,1],'r-', marker='.')
        #plt.plot(X,abs(U_t[k,:,1]-U_2t[k,:,1]),'k-')
        #plt.plot(X,U_t[k,:,1]/U_t[k,:,0],'r-')
    if(nb_var>2):
        1+1
        plt.plot(X,U_t[k,:,2],'g-')
        plt.plot(X,U_2t[k,:,2],'g-', marker='.')
        #plt.plot(X,abs(U_t[k,:,2]-U_2t[k,:,2]),'k-')
    
    #plt.plot(X,P_[k,:],'k-x')
    m=-.2; M=1.2
    
    #plt.ylim(3.8,4.8); plt.xlim(0.2,0.8)
    plt.ylim(m,M)
    #plt.xlim(-3.8, 3.8)
    plt.show()    
    #plt.savefig("save_"+str(T[k])+".png")    
    #plt.cla() 
 
    
if(False):
    with open('buckley.dat') as f:
        lines = f.readlines()
        
    X_B   = np.zeros(30000)
    U_B   = np.zeros(30000)
    for i in range(0,30000):
        X_B[i] = lines[i][0:16]
        U_B[i] = lines[i][17:32]
        
    plt.plot(X_B,U_B, 'r')
    
    plt.plot(X,U_t[k],'b', marker=".")
    #plt.plot(X,U_t2[k],'go')
    plt.ylim(-.2,1.2)
    #plt.show()     
    plt.savefig("endfig_save.png")    


    