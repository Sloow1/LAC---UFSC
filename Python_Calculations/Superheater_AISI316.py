# -*- coding: utf-8 -*-
"""
Created on Wed May 11 11:07:52 2022

@author: Eduardo Forster Beathalter

Código para dimensionamento de um trocador de calor pelo método e-NUT.
Geometria baseada num duto helicoidal de seção circular envolto por paredes cilíndricas tangentes à hélice.
"""

import numpy as np

def secant(f,x0,x1,TOL = 1E-3,NMAX = 100):
    """
        Aplica o método da secante para encontrar a raíz de uma função f(x).
        
        Parameters
        ----------
            f : function
                Função que se deseja obter a raíz
            x0 : float
                Estimativa inicial 1
            x1 : float
                Estimativa inicial 2
            TOL : float, optional
                Tolerância mínima exigida
            NMAX : int, optional
                Número máximo de iterações
        
        Returns
        -------
        int 
            Número de iterações final
        float
            Raíz aproximada da função
                
    """
    assert x0 != x1, "Os valores iniciais precisam ser diferentes"
    c = 0
    x2 = 0.0
    while(abs(x1 - x0)>=TOL):
        x2 = x1 - ((f(x1)*(x1-x0))/(f(x1)-f(x0)))
        x0 = x1
        x1 = x2
        c += 1
        assert c < NMAX, "O máximo de iterações foi alcançado"
    return c, x2

mDotAir = 0.003179 #kg/s
mDotSteam = 0.000476737 #kg/s
cp_air = 1106 #J/kg.K, ar a 831.6 K, interpolado
cp_steam = 1991 #J/kg.K, vapor a 523.15 K, interpolado
T_in_air = 600+273.15 #K
T_in_steam = 100 + 273.15 #K
T_out_steam = 400 + 273.15 #K
C_air = mDotAir*cp_air #W/K
C_steam = mDotSteam*cp_steam #W/K
T_out_air = T_in_air - C_steam/C_air*(T_out_steam-T_in_steam) #K

if(C_air>C_steam):
    Cmin = C_steam
    Cr = Cmin/C_air
elif(C_steam>C_air):
    Cmin = C_air
    Cr = Cmin/C_steam
else:
    C_min = C_air
    Cr = 1
Qmax = abs(Cmin*(T_in_air - T_in_steam))#W
Q = mDotSteam*cp_steam*(T_out_steam-T_in_steam)#W
epsilon = Q/Qmax #-

def e(NUT):
    #Relação e-NUT para trocador casco tubo com 1 passe no casco
    return 2*(1+Cr+(1+Cr**2)**0.5*(1+np.exp(-NUT*(1+Cr**2)**0.5))/(1-np.exp(-NUT*(1+Cr**2)**0.5)))**-1
def e1(NUT):
    #Função para encontrar-se a raíz
    return e(NUT)-epsilon
NUT = secant(e1, 1, 2, 1e-6)[1]
UA = NUT*Cmin #W/K

de = 10.29e-3 #m, diâmetro externo do duto
L = 1.24e-3 #m, espessura da parede do duto
di = de - 2*L #m, diâmetro interno do duto
a = 20e-3 #m, distância entre centros dos dutos, passo da hélice

rho_air = 0.4192 #kg/m^3, ar a 831.6 K, interpolado
rho_steam = 0.42198 #kg/m^3, vapor a 523.15 K, interpolado
A_air = de*(a-np.pi*de/4) #m^2
A_steam = np.pi*(di**2)/4 #m^2
V_air = mDotAir/rho_air/A_air #m/s
V_steam = mDotSteam/rho_steam/A_steam #m/2
visc_air = 379e-7 #Pa.s, ar a 831.6 K, interpolado
visc_steam = 178.7e-7 #Pa.s, vapor a 523.15 K, interpolado

if(a-de<de):
    Re_air = rho_air*V_air*(a-de)/visc_air
else:
    Re_air = rho_air*V_air*de/visc_air
print("Velocidade do ar: "+str(V_air)+" m/s")
print("Velocidade do vapor: "+str(V_steam)+" m/s")
Re_steam = rho_steam*V_steam*di/visc_steam
print("Reynolds do ar: "+str(Re_air))
print("Reynolds do vapor: "+str(Re_steam))

f_air = (0.79*np.log(Re_air)-1.64)**-2 #-
f_steam = (0.79*np.log(Re_steam)-1.64)**-2 #-
Pr_air = 0.713 #Prandtl ar a 831.6 K, interpolado
Pr_steam = 0.996 #Prandtl vapor a 523.15 K, interpolado
Nu_air = ((f_air/8)*(Re_air-1000)*Pr_air)/(1+(12.7*(f_air/8)**0.5)*(Pr_air**(2/3)-1))
k_air = 58.75e-3 #W/m.K, ar a 831.6 K, interpolado
Nu_steam = ((f_steam/8)*(Re_steam-1000)*Pr_steam)/(1+(12.7*(f_steam/8)**0.5)*(Pr_steam**(2/3)-1))
k_steam = 35.75e-3 #W/m.K, vapor a 523.15 K, interpolado
h_air = Nu_air*k_air/L #W/m^2.K
h_steam = Nu_steam*k_steam/di #W/m^2.K

k = 13.4 #W/m.K, para AISI 316
Rcond = L/k #m^2.K/W
Rconv_air = 1/h_air #m^2.K/W
Rconv_steam = 1/h_steam #m^2.K/W
Rtot = Rconv_air + Rconv_steam + Rcond #m^2.K/W
Utot = 1/Rtot #W/m^2.K

print("Coeficiente de transferência global: "+str(Utot)+" W/m^2.K")

Atot = UA/Utot #m^2
r = 52e-3 #m, raio da espira
ne = Atot/np.pi/de/(((2*np.pi*r)**2+a**2)) #Nr de espiras da serpentina
print("Temperatura de saída do vapor: "+str(T_out_steam-273.15)+" ºC")
print("Temperatura de saída do ar: "+str(T_out_air-273.15)+" ºC")
print("Número de espiras: "+str(ne))
