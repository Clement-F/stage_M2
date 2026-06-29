import numpy as np
import matplotlib.pyplot as plt


p = lambda U:  0.4*(U[2] - 0.5*(U[1]*U[1])/U[0])

with open('file_data.txt') as f:
    lines = f.readlines()
nb_var =int(lines[0][10:16]) 
order = int(lines[1][10:16])
nx = int(lines[3][5:11])
sm = int(lines[4][5:11])

T=[]
#for k in range(4,4+sm):
 #   T.append(float(lines[k][14:30]))


#   solDG   SolSub   file_sol    

with open('solDG.txt') as f:
    lines2 = f.readlines()   

with open('solSub.txt') as f:
    lines1 = f.readlines()   

U_t = np.zeros((sm,nx,nb_var)); P_ = np.zeros((sm,nx))
U_2t = np.zeros((sm,nx,nb_var)); P2_ = np.zeros((sm,nx))
X   = np.zeros(nx)
dec = 0


for k in range(0,sm):
    for i in range(0,(nx)): 
        X[i] = lines1[k*(nx+1) + i][0:10]
        U_t[k][i][0]   = lines1[k*(nx+1) +i][11:27]
        if(nb_var >1) : 
            U_t[k][i][1]   = lines1[k*(nx+1) +i][28:43]
        if(nb_var >2) :
            U_t[k][i][2]   = lines1[k*(nx+1) +i][44:58]
            
            
        U_2t[k][i][0]   = lines2[k*(nx+1) +i][11:27]
        print(max(abs(U_t[k,:,0]-U_2t[k,:,0])))
        if(nb_var >1) : 
            U_2t[k][i][1]   = lines2[k*(nx+1) +i][28:43]
            print(max(abs(U_t[k,:,1]-U_2t[k,:,1])))
        if(nb_var >2) :
            U_2t[k][i][2]   = lines2[k*(nx+1) +i][44:58]
            print(max(abs(U_t[k,:,2]-U_2t[k,:,2])))
            
        
        #P_[k,i] = p(U_t[k,i,:])
                  
        
    plt.plot(X,U_t[k,:,0],'b-')
    plt.plot(X,U_2t[k,:,0],'b-')
    plt.plot(X,abs(U_t[k,:,0]-U_2t[k,:,0]),'k-')
    if(nb_var>1):
        plt.plot(X,U_t[k,:,1],'r-')
        plt.plot(X,U_2t[k,:,1],'r-')
        plt.plot(X,abs(U_t[k,:,1]-U_2t[k,:,1]),'k-')
        #plt.plot(X,U_t[k,:,1]/U_t[k,:,0],'r')
    if(nb_var>2):
        plt.plot(X,U_t[k,:,2],'g')
        plt.plot(X,U_2t[k,:,2],'g')
        plt.plot(X,abs(U_t[k,:,2]-U_2t[k,:,2]),'k-')
    
    #plt.plot(X,P_[k,:],'k')
    
    plt.ylim(-1.2,1.2)
    plt.show()         


    