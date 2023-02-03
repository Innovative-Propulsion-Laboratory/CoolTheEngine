# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 20:47:01 2021

@author: julie
"""

from casthermo import *
import numpy as np
from matplotlib import pyplot as plt
import scipy as sc
"""
pas=0.015
epaisseur=0.001
hauteur=0.0018
largeur=0.0018
dx=0.00005
Hg=2600
lamb=388
Tg=3200
Hl=59500
Tl=158.5

pas=0.015
epaisseur=0.001
hauteur=0.0015
largeur=0.0012
dx=0.0001
Hg=8000
lamb=390
Tg=3000
Hl=120000
Tl=125
"""
def carto2D(pas,epaisseur,hauteur,largeur,dx,Hg,lamb,Tg,Hl,Tl,w,oui,leg):
    #def de la zone paroi
    npxp=round(pas/(2*dx)+1,0)
    npyp=round((epaisseur)/dx+1,0)
    #print(npxp,npyp)
    #def de la zone ailette
    npxa=round(((pas-largeur)/(2*dx))+1,0)
    npya=round((hauteur)/dx,0)
    #print(npxa,npya)
    #def du nombre de point à chercher
    nbp=npxp*npyp+npxa*npya
    #print(nbp)
    
    #def du type de config, orientation, symétrie et indice des points alentours
    listing=[]
    coord=[]
    #rangée de surface
    for i in range (0,int(npxp),1):
        point=[2,3]
        if i==0:
            sym=1
            n=-1
            e=len(listing)+1
            s=len(listing)+npxp
            o=-1
        elif i==(npxp-1):
            sym=2
            n=-1
            e=-1
            s=len(listing)+npxp
            o=len(listing)-1    
        else:
            sym=0
            n=-1
            e=len(listing)+1
            s=len(listing)+npxp
            o=len(listing)-1 
        u=i*dx
        v=0
        coord.append([u,v])
        card=[n,o,s,e]
        point.append(sym)
        point.append(card)
        listing.append(point)
    #rangée dans la paroi
    for h in range (0,int(npyp)-2,1):
        for i in range (0,int(npxp),1):
            point=[3,0]
            if i==0:
                sym=1
                o=-1
                n=len(listing)-npxp
                e=len(listing)+1
                s=len(listing)+npxp
            elif i==(npxp-1):
                sym=2
                o=len(listing)-1
                n=len(listing)-npxp
                e=-1
                s=len(listing)+npxp
            else:
                sym=0
                o=len(listing)-1
                n=len(listing)-npxp
                e=len(listing)+1
                s=len(listing)+npxp
            card=[n,o,s,e]
            point.append(sym)
            point.append(card)
            listing.append(point) 
            u=i*dx
            v=h*dx+dx
            coord.append([u,v])
    #changement de fluide
    abop=len(listing)-1
    for i in range (0,int(npxp),1):
        if i==0:
            point=[2,1]
            e=len(listing)+1
            o=-1
            n=len(listing)-npxp
            s=-1
        elif i<=(npxp-npxa-1):
            point=[2,1]
            e=len(listing)+1
            o=len(listing)-1
            n=len(listing)-npxp
            s=-1
        elif i==(npxp-npxa):
            point=[2,1]
            e=len(listing)+1
            o=len(listing)-1
            n=len(listing)-npxp
            s=len(listing)+npxa
        elif i==(npxp-1):
            point=[3,0]
            e=-1
            o=len(listing)-1
            n=len(listing)-npxp
            s=len(listing)+npxa
        else:
            point=[3,0]
            e=len(listing)+1
            o=len(listing)-1
            n=len(listing)-npxp
            s=len(listing)+npxa
        if i==0:
            sym=1
        elif i==(npxp-1):
            sym=2
        else:
            sym=0
        card=[n,o,s,e]
        point.append(sym)
        point.append(card)
        listing.append(point)
        u=i*dx
        v=(npyp-1)*dx
        coord.append([u,v])
    for h in range (0,int(npya)-1,1):
        for i in range (0,int(npxa),1):
            if i==0:
                point=[2,2]
                o=-1
                n=len(listing)-npxa
                e=len(listing)+1
                s=len(listing)+npxa
            elif i==(npxa-1):
                point=[3,0]
                o=len(listing)-1
                n=len(listing)-npxa
                e=-1
                s=len(listing)+npxa         
            else:
                point=[3,0]
                o=len(listing)-1
                n=len(listing)-npxa
                e=len(listing)+1
                s=len(listing)+npxa
            if i==(npxa-1):
                sym=2
            else:
                sym=0
            card=[n,o,s,e]
            point.append(sym)
            point.append(card)
            listing.append(point)
            u=i*dx+(npxp-npxa)*dx
            v=(npyp-1)*dx+(h+1)*dx
            coord.append([u,v])
    for i in range (0,int(npxa),1):
        if i==0:
            point=[4,0]
            o=-1
            n=len(listing)-npxa
            e=len(listing)+1
            s=-1
        elif i==(npxa-1):
            point=[5,0]
            o=len(listing)-1
            n=len(listing)-npxa
            e=-1
            s=-1
        else:
            point=[5,0]
            o=len(listing)-1
            n=len(listing)-npxa
            e=len(listing)+1
            s=-1
        if i==(npxa-1):
            sym=2
        else:
            sym=0
        card=[n,o,s,e]
        point.append(sym)
        point.append(card)
        listing.append(point)
        u=i*dx+(npxp-npxa)*dx
        v=(npya+npyp-1)*dx
        coord.append([u,v])
    #print(len(listing))
    #print(listing)
    
    ###résolution de la matrice inversible
    #shaping de la matrice inversible
    reso=np.zeros(shape=(int(nbp),int(nbp)))
    membre=np.zeros(shape=(int(nbp),1))
    for k in range (0,len(listing),1):
    #définition du fluide si présent en ce point
        if k<abop:
            h=Hg
            Tf=Tg
            inv=1
        else:
            h=Hl
            Tf=Tl
            inv=1
    #résolution des coefficients
        if listing[k][0]==1:
            a,b,c,d,x,plus=cas1(h,dx,lamb,Tf,inv)
        elif listing[k][0]==2:
            a,b,c,d,x,plus=cas2(h,dx,lamb,Tf,inv)
        elif listing[k][0]==3:
            a,b,c,d,x,plus=cas3(h,dx,lamb,Tf,inv)        
        elif listing[k][0]==4:
            a,b,c,d,x,plus=cas4(h,dx,lamb,Tf,inv)    
        elif listing[k][0]==5:
            a,b,c,d,x,plus=cas5(h,dx,lamb,Tf,inv)    
    #résolution de l'orientation  
        if listing[k][1]==0:
            coef1=a
            coef2=b
            coef3=c
            coef4=d
        elif listing[k][1]==1:
            coef1=b
            coef2=c
            coef3=d
            coef4=a
        elif listing[k][1]==2:
            coef1=c
            coef2=d
            coef3=a
            coef4=b
        elif listing[k][1]==3:
            coef1=d
            coef2=a
            coef3=b
            coef4=c
        
    #résolution symétrie
        if listing[k][2]==1:
            coef4=coef2+coef4
        elif listing[k][2]==2:
            coef2=coef4+coef2
    
    #déduction du placement des coef par rapport aux points accolant
        insert=[]
        pos=1
        for z in listing[k][3]:
            if z>=0 and z<=(nbp-1) :
                if pos==1:
                    inn=[coef1,z]
                    insert.append(inn)
                elif pos==2:
                    inn=[coef2,z]
                    insert.append(inn)
                elif pos==3:
                    inn=[coef3,z]
                    insert.append(inn)
                elif pos==4:
                    inn=[coef4,z]
                    insert.append(inn)
                else:
                    print("error placing")
            pos=pos+1
        
    #introduction des coefs dans la matrice
        for values in insert :
            implacement=int(values[1])
            reso[k][implacement]=values[0]
            reso[k][k]=x
            membre[k]=plus
    #print(reso)
    #print(membre)
    reso_inv=np.linalg.inv(reso)
    #reso_inv=sc.linalg.inv(reso)
    #print(np.dot(reso_inv,reso))
    T=np.dot(reso_inv,membre)
    #print(T)
    
    minimum=min(T)
    maximum=max(T)
    
    abcisse=[]
    ordonee=[]
    temperature=[]
    bis=1
    paroigas=0
    for m in coord:
        abcisse.append(m[0])
        ordonee.append(-m[1])
    for t in T:
        temperature.append(t[0])
        if bis<=npxp:
            paroigas=paroigas+t[0]
        bis=bis+1
    if oui==1:
        print("█ Température moyenne à la paroi gas =",round(temperature[0],3),"                            █")  
        print("█ Température moyenne paroi coolant  =",round(temperature[abop+1],3),"                            █")
        print("█ Température maximum à la paroi     =",round(max(temperature),5),"                            █")  
        if leg==1:
            a1=0.004
            a2=0.0005
            a3=0.00025
            a4=-0.002
        elif leg==2:
            a1=0.002
            a2=0.00025
            a3=0.00025
            a4=-0.0015
        plt.figure(dpi=200)
        p=plt.scatter(abcisse, ordonee ,c=temperature ,marker='s',s=w, cmap='flag')#rainbow#prism#flag
        plt.text(a1, a2, 'Hot gases', horizontalalignment = 'center', verticalalignment = 'center')
        plt.text(a3, a4, 'Coolant', horizontalalignment = 'center', verticalalignment = 'center')

        plt.title("Distribution 2D des températures",fontsize=15)
        plt.axis("equal")
        plt.colorbar(p,shrink=0.4, aspect=15)
        plt.show()
        
        plt.figure(dpi=200)
        p=plt.scatter(abcisse, ordonee ,c=temperature ,marker='s',s=w, cmap='rainbow')#rainbow#prism#flag
        plt.text(a1, a2, 'Hot gases', horizontalalignment = 'center', verticalalignment = 'center')
        plt.text(a3, a4, 'Coolant', horizontalalignment = 'center', verticalalignment = 'center')

        plt.title("Distribution 2D des températures",fontsize=15)
        plt.axis("equal")
        plt.colorbar(p,shrink=0.4, aspect=15)
        plt.show()
    
    t3d=[]
    for i in range (0,int(npxp),1):
        t3d.append(temperature[i])
    
    return t3d
        
        

#r=carto2D(pas,epaisseur,hauteur,largeur,dx,Hg,lamb,Tg,Hl,Tl,20)

    
            









        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    

