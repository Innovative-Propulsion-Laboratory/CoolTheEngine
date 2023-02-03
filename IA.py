# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 23:47:42 2022

@author: julie
"""
from random import *
import numpy as np
from math import *
from sympy import Symbol, nsolve, Eq
import sympy as mp

def AlgoDE(D,NP,F,CR,Gmax,xmin,xmax,Gdeb,HHV,Xqual,We,Kp,Xeq,condcoolant,Dhy,ycanauxre,epaiss_chemise,htre,nbc,c,Lambda_tc,Tcoolant,hg,Tg,qr):
    def function(flux1,Gdeb,HHV,Xqual,We,Kp,Xeq,condcoolant,Dhy,ycanauxre,epaiss_chemise,htre,nbc,c,Lambda_tc,Tcoolant,hg,Tg,qr):
        Bo = 1000000*flux1/(Gdeb*HHV)
        #print(Xqual,Xeq)
        if Xqual<0.6:
            Nu = 12.46*(Bo**0.544)*(We**0.035)*(Kp**0.614)*(Xeq**0.031)
            #print("A",Nu)
        else:
            Nu = 0.00136*(Bo**(-1.442))*(We**0.074)
            #print("B",flux1,Gdeb,HHV)
        #print(Nu)
        hl=Nu*(condcoolant/Dhy)
        #print(" ")
        #print(Nu,condcoolant[i],Dhy)
        #hlnormal.append(hl)
        D=2*(ycanauxre-epaiss_chemise)
        d=(pi*(D+htre+epaiss_chemise)-nbc*c)/nbc
        #print(" ")
        #print(d,hl,Lambda_tc)
        #print("l'erreur est ",V,"ou",Dhy,"ou",rho[i],"ou",visccoolant[i],Tcoolant[i])
        m=((2*hl)/(d*Lambda_tc))**0.5
        #print(m,htre[i])
        #print(hl,nbc,c,pi,D,Lambda_tc,m,htre[i])
        #print(hl,nbc,c,D,m,htre[i])
        hl_cor=hl*((nbc*c)/(pi*D))+nbc*((2*hl*Lambda_tc*(((pi*D)/nbc)-c))**0.5)*((tanh(m*htre))/(pi*D))
        #hlcor.append(hl_cor)
        #Résolution températures paroi
        """
        hl=hl_cor
        Tl=Tcoolant #K
        e=epaiss_chemise
        L=Lambda_tc
        mp.dps = 150
        cx1 = Symbol('cx1')
        cx2 = Symbol('cx2')
        #print(hg,Tg,L,e,qr,hl,Tl)
        f1 = hg*(Tg-cx1)-(L/e)*(cx1-cx2)+qr
        f2 = hl*(cx2-Tl)-(L/e)*(cx1-cx2)
        x_,y_=nsolve((f1,f2), (cx1, cx2),(900,700))
        """
        y_=((Lambda_tc/epaiss_chemise)*(hg*Tg+hl_cor*Tcoolant)+hl_cor*hg*Tcoolant)/(hg*hl_cor+(Lambda_tc/epaiss_chemise)*(hg + hl_cor))
        x_=((Lambda_tc/epaiss_chemise)*(hg*Tg+hl_cor*Tcoolant)+hl_cor*hg*Tg)/(hg*hl_cor+(Lambda_tc/epaiss_chemise)*(hg + hl_cor))
                
        #print(x_,y_)
        #inwall_temperature.append(x_)
        #outwall_temperature.append(y_)
        
        #calcul du flux
        flux2=hl*(y_-Tcoolant)*0.000001
        return (abs(flux1-flux2)/flux2)       
    #DE Algorithm    
    G=0     #Generations
    P=np.zeros((D,NP))
    Y=np.zeros((Gmax,NP))
    NPbis=[]
    Gfollow=[0]
    for i in range (0,NP,1):
        for j in range(0,D,1):
            x_de=round(uniform(xmin[j],xmax[j]),6)
            P[j,i]=x_de
            #NPbis[j,i]=i
        y_de=function(P[0,i],Gdeb,HHV,Xqual,We,Kp,Xeq,condcoolant,Dhy,ycanauxre,epaiss_chemise,htre,nbc,c,Lambda_tc,Tcoolant,hg,Tg,qr)
        Y[G,i]=y_de
        #print(sec,"Generation",G,"Indiv.",i,"Rate",Y[G,i],P[0,i],P[1,i],P[2,i],P[3,i],P[4,i],P[5,i],P[6,i],P[7,i],P[8,i],P[9,i],P[10,i],P[11,i],P[12,i],P[13,i])         
    G=G+1
    for jj in range (0,D,1):
        LL=[]
        for ii in range (0,NP,1):
            LL.append(ii)
        NPbis.append(LL)
    #print(NPbis)
    #best=Y.index(min(Y))
    #Xprim[0,:]=P
    y1cost=1
    while G<Gmax and y1cost>0.001:
        Pmove=np.zeros((D,NP))
        Ymove=np.zeros((1,NP))
        for s in range (0,NP,1):
            x_i=np.zeros((D,1))
            for d in range(0,D,1):
                x_i[d,0]=P[d,s]
                NPbis[d].remove(s)
            #print(NPbis) 
            #Mutation
            v=[]
            for d in range(0,D,1):
                cho=[]
                for r in range (0,3,1):
                    r_all=choice(NPbis[d])
                    NPbis[d].remove(r_all)
                    cho.append(r_all)
                    #print(cho)
                P_r1=P[d,cho[0]]
                P_r2=P[d,cho[1]]
                P_r3=P[d,cho[2]]
                v_=P_r1+F*(P_r2-P_r3)
                v_=round(v_,6)
                v.append(v_)
                for r in range (0,3,1):
                    NPbis[d].append(cho[r])
                NPbis[d].append(s)     
            #Crossover
            u=np.zeros((D,1))
            for d in range(0,D,1):
                CR_bis=uniform(0,1)
                if CR_bis<=CR and v[d]<=xmax[d] and v[d]>=xmin[d]:
                    u[d,0]=v[d]
                else:
                    u[d,0]=x_i[d,0]     
            #Comparaison
            y1cost=function(u[0,0],Gdeb,HHV,Xqual,We,Kp,Xeq,condcoolant,Dhy,ycanauxre,epaiss_chemise,htre,nbc,c,Lambda_tc,Tcoolant,hg,Tg,qr)           
            if y1cost<Y[G-1,s]:
                Pmove[:,s]=u[:,0]
                Ymove[0,s]=y1cost
            else:
                Pmove[:,s]=x_i[:,0]
                Ymove[0,s]=Y[G-1,s]
            #print(sec,"Generation",G,"Indiv.",s,"Rate",Ymove[0,s],Pmove[0,s],Pmove[1,s],Pmove[2,s],Pmove[3,s],Pmove[4,s],Pmove[5,s],Pmove[6,s],Pmove[7,s],Pmove[8,s],Pmove[9,s],Pmove[10,s],Pmove[11,s],Pmove[12,s],Pmove[13,s])    
        P=Pmove
        Y[G,:]=Ymove[0,:]
        #Gfollow.append(G+1)
        #print(G)
        #Xprim[G+1,:]=P[:]
        G=G+1
        Gfollow.append(G)
    print(P,min(Y[-1,:]))
    return P[0,0]