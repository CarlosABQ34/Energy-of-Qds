import math as ma
import numpy as np
import scipy.special as sc
import csv 
import time
from multiprocessing import Pool
from os import cpu_count
def Calculos(x: float):
    #Numeric values
    inicio = time.time() #Function for counting
    e = 37.9407 #charge mev^(1/2)*nm^(1/2)
    e_2=1439.5 #charge to square
    #CdTe
    mo=5.674*pow(10,-9)
    me=0.096*mo
    mh=0.4*mo
    k=10.4
    gap=1610
    h = 6.58211899*pow(10,-4)  # hbarmeV*ns.
    #parameters
    '''
    CdTe
    me=0.096*mo    mh=0.4*mo    k=10.4    gap=1610
    CdSe
    me=0.13*mo    mh=0.4*mo    k=10.6    gap=1750
    ZnTe
    me=0.15*mo    mh=0.2*mo    k=9.7    gap=2200
    '''
    #radius confinament(nm)
    asp=x
    r=7
    volumen=4*ma.pi*pow(r,3)/3 #(nm**3)
    '''
    We take values such that lh and lz are equal for both the hole and the electron, with subscripts _1 denoting the hole and _2 denoting the electron.
    '''
    def lzlh(r, asp):
        lz =r*pow(asp,2/3) 
        lh = r/pow(asp,1/3)
        return lz, lh
    lz, lh =lzlh(r, asp)
    #Frecuencus for electron and hole
    wze=h/(2*me*pow(lz,2))
    wle=h/(2*me*pow(lh,2))
    wzh=h/(2*mh*pow(lz,2))
    wlh=h/(2*mh*pow(lh,2))
    xe=pow(lh,2)/pow(lz,2)
    #Energy oscillator
    def Mango(n1p,m1p,s1p,n2p,m2p,s2p,n1,m1,s1,n2,m2,s2):
        if (n1==n1p and n2==n2p and m1==m1p and m2==m2p and s1==s1p and s2==s2p):
            E=h*wle*(n1+0.5)+h*wlh*(n2+0.5)+h*wle*(m1+0.5)+h*wlh*(m2+0.5)+h*wze*(s1+0.5)+h*wzh*(s2+0.5)+gap
        else:
            E=0
        return E
    #Interaction function 
    def Vco(n1p,m1p,s1p,n2p,m2p,s2p,n1,m1,s1,n2,m2,s2): 
        Vcoul=0
        Rr=m1+m2-(n1+n2)
        Rl=m1p+m2p-(n1p+n2p) 
        Rs=s1+s1p+s2+s2p
        Unid=e_2/(2*ma.pi*k*lh)
        raiz=1/(ma.sqrt(sc.factorial(n1)*sc.factorial(n1p)*sc.factorial(n2)*sc.factorial(n2p)*sc.factorial(m1)*sc.factorial(m1p)*sc.factorial(m2)*sc.factorial(m2p)*sc.factorial(s1)*sc.factorial(s1p)*sc.factorial(s2)*sc.factorial(s2p)))
        if (Rl==Rr and Rs%2==0):
            for p1 in range(0,min(n1,n1p)+1):
                for p2 in range(0,min(m1,m1p)+1):
                    for p3 in range(0,min(s1,s1p)+1):
                        for p4 in range(0,min(n2,n2p)+1):
                            for p5 in range(0,min(m2,m2p)+1):
                                for p6 in range(0,min(s2,s2p)+1):
                                    u=m1p+m2p+n1+n2-(p1+p2+p4+p5)
                                    v=s1+s1p+s2+s2p-2*(p3+p6)
                                    a1=1+u
                                    a2=(1+v)*0.5
                                    b=(1+2*u+v)*0.5
                                    c=(3+2*u+v)*0.5             
                                    if (asp>=1):
                                        f1=pow(xe,u+0.5)*sc.hyp2f1(a1,b,c,1-xe)
                                    else:
                                        f1= pow(xe,-v*0.5)*sc.hyp2f1(a2,b,c,1-(1/xe)) 
                                    Vcoul+=Unid*sc.factorial(p1)*sc.factorial(p2)*sc.factorial(p3)*sc.factorial(p4)*sc.factorial(p5)*sc.factorial(p6)* sc.comb(n1,p1)*sc.comb(n1p,p1)*sc.comb(m1,p2)*sc.comb(m1p,p2)*sc.comb(s1,p3)*sc.comb(s1p,p3)*sc.comb(n2,p4)*sc.comb(n2p,p4)*sc.comb(m2,p5)*sc.comb(m2p,p5)*sc.comb(s2,p6)*sc.comb(s2p,p6)*raiz*f1*pow(-1,u+v*0.5+n2p+m2p+s2p+n2+m2+s2)*pow(0.5,u)*(sc.gamma(a1)*sc.gamma(a2)*sc.gamma(b))/sc.gamma(c)  
        else:
            Vcoul=0
        return Vcoul  
    #Quantum numbers lecture
    rows = np.loadtxt('tablas2.csv', delimiter=',', dtype=int)
    m = len(rows)
    mij = np.zeros((m, m))
    Eij = np.zeros((m, m))
    primer=np.zeros((m, m))
    #Matrix calculation
    for i in range (m):
        for j in range (m):
            Eij[i,j]=Mango(*rows[i], *rows[j])#Oscillar matrix calculation
    for i in range(m):
        for j in range(m):
            mij[i,j]=Vco(*rows[i], *rows[j])#Ic matrix 
    for i in range(m):
        primer[i, i] = Vco(*rows[i], *rows[i])  # Asignar el valor de la función en la diagonal
       #Sin interaccion        
    print("Matriz sin inte")
    print(Eij)
    #interaccion
    print("matriz pertubacion")
    print(primer)
    print("matriz de IC")
    print(mij)
    #Saving eigenvalues without interaction
    auValores, auVectores = np.linalg.eigh(Eij)
    auv=auValores[:30]
    print(auv)
    file_name=f"Orauto_{asp:1.2f}lh{lh:1.2f}lz{lz:1.2f}r{r}.txt"
    with open (file_name, "w") as temp_file:
        for item in auv:
            temp_file.write("%s\n" % item)
    file=open(file_name,"r")
    #Saving eigenvalues in 1st pertubation matrix
    G=Eij-primer
    valp, vectp = np.linalg.eigh(G)
    valp5=valp[:30]
    np.set_printoptions(precision=6, suppress=True)
    file_name3=f"Valorp{asp:1.2f}lh{lh:1.2f}lz{lz:1.2f}r{r}.txt"
    with open (file_name3, "w") as temp_file:
        for item in valp5:
            temp_file.write("%s\n" % item)
    file=open(file_name3,"r")
    #Total Hamiltonian
    Hami=Eij-mij
    print("Hamiltoniano Total")
    print(Hami) 
    autva, autve = np.linalg.eigh(Hami)
    #Saving eigenvectors with interaction
    autva5=autva[:30]
    np.set_printoptions(precision=6, suppress=True)
    #Saving eigenvalues without interaction
    file_name3=f"Autova{asp:1.2f}lh{lh:1.2f}lz{lz:1.2f}r{r}.txt"
    with open (file_name3, "w") as temp_file:
        for item in autva5:
            temp_file.write("%s\n" % item)
    file=open(file_name3,"r")
    np.set_printoptions(precision=6, suppress=True) 
    print ("Los primeros 50 autovalore son")
    print(autva5)
    fin = time.time()
    tiempo_transcurrido = fin - inicio
    print("Tiempo transcurrido:", tiempo_transcurrido, "segundos")
#Main function
if __name__ == "__main__":
    #aspect ratios (0,5-2)
    x = []
    for k in range(21):
        if k < 10:
            x.append(k / 20 + 0.5)
        elif k == 10:
            x.append(1.0)
        else:
            value = round(float(x[k - 1]) + 0.1, 1)
            x.append(value if value <= 2 else '2')
    print(x)
    with Pool(cpu_count()) as p:
        p.map(Calculos,x)