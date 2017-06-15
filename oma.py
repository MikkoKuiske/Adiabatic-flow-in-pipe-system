# -*- coding: utf-8 -*-
"""
Created on Wed Jun 07 14:11:58 2017

@author: pirttinm
"""
from __future__ import division
import numpy as np
from scipy.optimize import fsolve
from CoolProp.CoolProp import PropsSI

def main():
    
    ### input
    m1 = 8.63e-4
    T1 = 20.3+273.15
    p1 = 1e5
    d = 3e-3
    l = 50*d
    ### static constants
    R = 287
    gamma = 1.4
    roughness = 0.015e-3
    
    
    ### constants
    rho = PropsSI("D","T",T1,"P",p1,"Air")
    mu = PropsSI("V","T",T1,"P",p1,"Air")
    nu = mu/rho
    f_guess = 0.005
    Ma2_guess = 0.5
    
    ### calculation
    V1 = m1/(np.pi*d**2/4*rho)
    Re = V1*d/nu
    
    
    # friction factor, colebrook-white
    f_solver = lambda f_solved: 1/np.sqrt(f_solved)+2*np.log(roughness/(3.7*d)+2.51/(Re*np.sqrt(f_solved)))
    
    ksi = fsolve(f_solver,f_guess)[0]*4    
    Ma1 = V1/np.sqrt(gamma*R*T1)    
    
    ma_solver = lambda Ma_solved: ksi/d*l-1/gamma*(1/Ma1**2-1/Ma_solved**2)-(gamma-1)/(2*gamma)*\
                np.log(Ma1**2/Ma_solved**2*(1+0.5*(gamma-1)*Ma_solved**2)/((1+0.5*(gamma-1)*Ma1**2)))
            
    Ma2 = fsolve(ma_solver,Ma2_guess)[0]

    #  
    #L_s = d/(ksi*gamma)*((gamma+1)/2*np.log((gamma+1)/2*Ma1**2/(1+(gamma-1)/2*Ma1**2))+1/Ma1**2-1)
    
    
    V2 = np.sqrt(gamma*R*T1)*Ma2
    T2 = T1*(1/(gamma-1)+Ma1**2/2)/(1/(gamma-1)+Ma2**2/2)
    p2 = p1*(T1/T2)**(gamma/(gamma-1))
    
    ### results
    print "d =",d,"m","massflow",m1 
    print "ksi".ljust(6), round(ksi,4)
    print "p1".ljust(6),round(p1), "Pa"
    print "p2".ljust(6),round(p2), "Pa"
    print "T1".ljust(6),round(T1,2), "K"
    print "T2".ljust(6),round(T2,2), "K"
    print "V1".ljust(6),round(V1,2), "m/s"
    print "V2".ljust(6),round(V2,2), "m/s"
    print "Ma1".ljust(6),round(Ma1,3)
    print "Ma2".ljust(6),round(Ma2,3)
    #print "L_s".ljust(6), L_s
    
    #print V1*np.pi*d**2/4*60000
    #print V2*np.pi*d**2/4*60000
    
    
if __name__ == "__main__":
    main()