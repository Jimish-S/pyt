# -*- coding: utf-8 -*-
"""
Created on Tue Mar 01 08:45:56 2016

@author: USER"""
""" An elementary gas phase reaction is carried out in packed bed reactor. There is heat exchanging jacket surrounding the instrument. Find concentration,
conversion, temperatureand reduced pressure profile along the length of catalyst packing.

NOTE: I haven't considered catalyst degradation. For time being, I have confined
myself to simple reaction of 2A <=> C . All are gases. In later stage I will 
extend the same principle to reactions of ETO with catalyst degradation and then
optimise the packing for given output. In final stage I will generalise this 
programme for any gaseous conversions by user input of reactant and initial
conditions."""

#importing the required packages
import scipy as py
import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import odeint 

#Data
Ta= 500.0 #K
delH = -40000.0 #J/kg
Cpa = 40.0 #J/gmolK
FaO = 5.0 #gmol/min

#Initial Condition
f_int = [0.0,450.0,1.0] # x,T,y where x= conversion,T = Temperature, y = Reduced pressure
w = range(0,21,1)

#Differential equations
def PBR(f,w):
    x,T,y = f
    k= 0.5*py.exp((41800.0/8.314)*(1.0/450.0-1.0/T)) #Rate constant 
    Ca = 0.271*(1.0-x)*(450.0/T)/(1.0-0.5*x)*y #Concentration of A
    Cc = 0.271*0.5*x*(450.0/T)/(1.0-0.5*x)*y #concentartion of C
    Kc = 25000.0*py.exp(delH/8.314*(1.0/450.0-1.0/T)) # equilibrium constant
    rA = -k*(Ca**2-Cc/Kc) #rate expression
    
    dxdw = -rA/FaO #design equation for PFR
    dTdw = (0.8*(Ta-T)+rA*delH)/(Cpa*FaO) #Heat balance
    dydw = -0.015*(1-0.5*x)*(T/450.0)/(2.0*y) #Pressure 
    return np.array([dxdw,dTdw,dydw],dtype = float)

#Solve differential equations
fs = odeint(PBR,f_int,w)
print fs
#plot solution
xs = fs[:,0]
Ts = fs[:,1]
ys = fs[:,2]
Ca = 0.271* (1-xs)*(450.0/Ts)/(1.0-0.5*xs)*ys
Cc = 0.271*0.5*xs*(450.0/Ts)/(1.0-0.5*xs)*ys

plt.plot(w,xs,'-',w,Ts/1000.0,'-.',w,ys,'.',w,Ca,'x',w,Cc,'+') #Temp is divided by 1000 to meet the scale
plt.title('Reaction conversion, temperature and concentrations')
plt.xlabel('w, catalyst weight, kg')
plt.legend(['x','T in K', 'y','Ca','Cc'],'best')
plt.show()
