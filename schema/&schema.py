import numpy as np
import matplotlib.pyplot as plt


p = lambda U:  0.4*(U[2] - 0.5*(U[1])**2 /U[0] )
    
with open('file_data.txt') as f:
    lines = f.readlines()
    
nb_var = int(lines[0][10:16])
orderx = int(lines[1][10:16])
ordert = int(lines[2][10:16])
nb_cell =int(lines[3][10:16]) ; nb_sub= int(lines[4][13:19])
nx = nb_cell *nb_sub*20

sm = int(lines[5][5:11])

T=[]
for k in range(6,6+sm):
   T.append(float(lines[k][14:30]))


#   solDG   SolSub   file_sol    

with open('solPol.txt') as f:
    lines1 = f.readlines()      
with open('solSub'+'ex'+'.txt') as f:
    lines2 = f.readlines()   
with open('solSub.txt') as f:
    lines3 = f.readlines()      


U_t = np.zeros((sm,nx,nb_var));  
U_2 = np.zeros((sm,nx,nb_var))   
U_ex = np.zeros((sm,nx,nb_var));  
X   = np.zeros(nx)

X_cell = np.zeros(nb_cell+1);       U_cell=  np.zeros((sm,nb_cell,nb_var))
X_midcell = np.zeros(nb_cell)
dec = 0

X_cell[-1]=1.


for k in range(0,sm,):
    for i in range(0,nx): 
        X[i] = lines1[k*(nx+1) + i][0:10]
        U_t[k][i][0]         = lines1[k*(nx+1) +i][10:27]
        U_2[k][i][0]         = lines3[k*(nx+1) +i][10:27]
        if(nb_var >1) : 
            U_t[k][i][1]   = lines1[k*(nx+1) +i][28:43]
        if(nb_var >2) :
            U_t[k][i][2]   = lines1[k*(nx+1) +i][44:58]
            
        if((i)%(nb_sub*20)==0) :
            X_cell[int((i)/(nb_sub*20))] =lines1[k*(nx+1) + i][0:10]
            U_cell[k][int((i)/(nb_sub*20))-1][0] = sum(U_t[k,i-nb_sub:i,0])/nb_sub
            X_midcell[int((i)/(nb_sub*20))] = (X_cell[int(i/(nb_sub*20))]+X_cell[int(i/(nb_sub*20))+1])/2
            
        # if(i%(orderx+1)==0) : 
        #     X_cell[int(i/(orderx+1))] = lines1[k*(nx+1) + i][0:10]
        #     U_cell[k][int(i/(orderx+1))-1][0] = sum(U_t[k,i-(orderx+1):i,0])/(orderx+1)
        #     X_midcell[int(i/(orderx+1))] = (X_cell[int(i/(orderx+1))]+X_cell[int(i/(orderx+1))+1])/2
            
    
        U_ex[k][i][0]   = lines2[k*(nx+1) +i][11:27]
        # if(nb_var >1) : 
        #      U_ex[k][i][1]   = lines2[k*(nx+1) +i][28:43]
        # if(nb_var >2) :
        #     U_ex[k][i][2]   = lines2[k*(nx+1) +i][44:58]
            
        
        #P_[k,i] = p(U_t[k,i,:])
    U_cell[k][-1][0] = sum(U_t[k,-(orderx+1)-1:-1,0])/(orderx+1)
           
    for i in range(0,nb_cell+1): 
        1+1
        plt.plot([X_cell[i],X_cell[i]],[-2,2],linestyle='--', color='gray')
        
    for i in range(0,nb_cell): 
        J=i*nb_sub*20
        plt.plot(X_cell[i],U_t[k,J,0],"b",marker="x")
        if(i>0): plt.plot(X_cell[i],U_t[k,J-1,0],"b",marker="x")        
        # plt.plot(X[J:J+orderx+1],U_t[k,J:J+orderx+1,0],'b-', marker='.')
    
    plt.plot(X_cell[-1],U_t[k,-1,0],"b",marker="x")        
    
    #plt.plot(X_midcell,U_cell[k,:,0], marker='.')
    #plt.plot(X,U_t[0,:,0],'g-')
    plt.plot(X,U_t[k,:,0],'b-',label="polynome")
    plt.plot(X,U_2[k,:,0],'r-',label="sous valeurs finis")
    plt.plot(X,U_ex[k,:,0],'g-',label="exacte")
    #plt.plot(X,abs(U_t[k,:,0]-U_ex[k,:,0]),'k-')
    #print(max(abs(U_t[k,:,0]-U_ex[k,:,0])))
    #print(sum(U_t[k,:,0]))
    #print(min(P_[k,:]))
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
    
    #plt.plot(X,P_[k,:],'k-.')
    #m=-.2; M=1.2
    
    #plt.ylim(3.8,4.8); plt.xlim(0.2,0.8)
   # plt.ylim(m,M);plt.xlim(-0.,1.)
    plt.ylim(-1.2, 1.2)
    #plt.title("plot at time "+str([T[k]]))
    plt.legend(loc='upper right')
    plt.show()    
    #plt.savefig("save_"+str(T[k])+".png")    
    #plt.cla() 
    
# plt.plot(X,U_t[k,:,0],'b', marker='.')
# plt.ylim(3.,4.8);plt.xlim(0.5,2.3)
# plt.savefig("save_"+str(T[k])+"zoom.png")    
    



    