#This code finds the angle between two normals of an ellipse
#and plots it.
import numpy as np
import matplotlib.pyplot as plt
import math

def ques(theta):
    A=np.array([3*np.cos(theta),math.sqrt(3)*np.sin(theta)])
    B=np.array([-3*np.sin(theta),math.sqrt(3)*np.cos(theta)])
    O=np.array([0,0])
    V1=np.array([1,0])
    V2=np.array([0,3])
    V=np.vstack((V1,V2))
    nA=np.matmul(V,A)              #direction vector of normal at A
    nB=np.matmul(V,B)              #direction vector of normal at B
    dotnA_nB=np.matmul(nA.T,nB)    #inner product of nA and nB
    Dr_2=np.square(np.linalg.norm(nA)*np.linalg.norm(nB))-np.square(dotnA_nB)
    Dr=math.sqrt(Dr_2)
    cot_beta=abs(float(dotnA_nB)/Dr)
    beta =math.atan(1/cot_beta)
    ans_nr=2*cot_beta
    ans_dr=-np.matmul(A.T,B)/float(3)
    ans=ans_nr/ans_dr
    print("Values of 2cotβ, sin2Θ and 2cotβ/sin2Θ are")
    print(ans_nr)
    print(ans_dr)  
    print (ans)
    a=3
    b=math.sqrt(3)
    t=np.linspace(0,2*math.pi,100)
    plt.plot(a*np.cos(t),b*np.sin(t))
    plt.grid()
    axis=plt.gca()
    plt.plot(A[0],A[1],'o')
    plt.plot(B[0],B[1],'o')
    plt.plot(O[0],O[1],'o')
    plt.text(A[0]*(1+0.1),A[1]*(1+0.1),'A')
    plt.text(B[0]*(1+0.1),B[1]*(1+0.1),'B')
    plt.text(0.1,0,'O')
    omat=np.array([[0,1],[-1,0]])
    n1=np.matmul(omat,nA)             #direction vector of tangent at A
    n2=np.matmul(omat,nB)             #direction vector of tangent at A
    n1_n2=np.vstack((n1,n2))
    p=np.zeros(2)
    p[0]=np.matmul(n1,A)
    p[1]=np.matmul(n2,B)
    X=np.matmul(np.linalg.inv(n1_n2),p)
    plt.plot(X[0],X[1],'o')
    plt.text(X[0]+0.1,X[1],'N')
    len=10
    lam_1=np.linspace(0,1,len)
    x_AX=np.zeros((2,len))
    x_BX=np.zeros((2,len))
    for i in range(len):
        temp1=A+lam_1[i]*(X-A)
        x_AX[:,i]=temp1.T
        temp2=B+lam_1[i]*(X-B)
        x_BX[:,i]=temp2.T
    plt.plot(x_AX[0,:],x_AX[1,:])
    plt.plot(x_BX[0,:],x_BX[1,:])
    plt.title("Normals to the ellipse at A and B intersects at N") 
    plt.axis([-4, 4, -2.5, 2.5])
    #polar graph plottation
    f, ax =plt.subplots(2, subplot_kw=dict(projection='polar'))
    ax[0].scatter(beta, 1)
    ax[0].set_title("β's                  graph");
    ax[1].scatter(theta, 1)
    ax[1].set_title("Θ's                  graph");
    f.subplots_adjust(hspace=0.5)
    #scatter point plottation
    f, ax1 =plt.subplots()
    ax1.set_title("Graph of 2cotβ vs sin2Θ")
    ax1.scatter(ans_nr, ans_dr)
    ax1.scatter(0,0)
    ax1.axis([-1.5, 1.5, -1.5, 1.5])
    
ques(np.pi/4) #put different values of theta
plt.legend(loc='best') 
plt.xlabel("x-axis")
plt.ylabel("y-axis")
plt.tick_params(axis='both', which='major',
labelsize=14)
plt.grid()
plt.show()
