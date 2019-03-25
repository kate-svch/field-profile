#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 12:23:46 2018

@author: kate-svch
"""
from math import*
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import time

import matplotlib as mpl
  # инициализация констант
epsilon0 = 8.85*10**(-12);  zLimit = 2000;
H=0; hm=9000; hu = 11000;    # hm and hu - heights of levels relative to LPCR!!!
#ql=11; qm=-60; qu=40;
sigma_l=700; sigma_m=1000; sigma_u=1000;
Rl=4000; Rm=12000; Ru=12000;
MountainHeight=2200;  

ro_l0=0; ro_m0=(-1.1)*10**(-9); 
#ro_u0=0.7*10**(-9); 
ro_u0 = 1.1*10**(-9);
ro_of_layers=[ ro_l0, ro_m0, ro_u0] 
ro_of_layers_Mirror=[ -ro_l0, -ro_m0, -ro_u0];    
sigma_of_layers=[ sigma_l, sigma_m, sigma_u]   
zCentre_of_layers=[H, H+hm, H+hu] 
zCentre_of_layers_Mirror=[-H, -(H+hm), -(H+hu)];
R_of_layers=[Rl, Rm, Ru]
 
 # задаём распределение плотности зарядов
# модуль выражения = 2*(выр.)-1
# это функция, задающая зависимость плотности заряда слоя номер j от координаты z
def Density_of_charge(j_layer, z, ro_of_layers):
 return ro_of_layers[j_layer] * exp( - ( ( z-zCentre_of_layers[j_layer] )/ sigma_of_layers[j_layer] )**8 )

def Density_of_charge_Mirror(j_layer, z, ro_of_layers_Mirror):
 return ro_of_layers_Mirror[j_layer] * exp( - ( ( z-zCentre_of_layers_Mirror[j_layer] )/ sigma_of_layers[j_layer] )**8 )

def zsignum(z, zInside):
    return 2*((-z+zInside)>0)-1;

def Elementary_field_function_z(j_layer, z, zInside, ro_of_layers):
 return Density_of_charge(j_layer, zInside, ro_of_layers) *(zsignum(z,zInside)+ (z-zInside)/sqrt( (R_of_layers[j_layer])**2 + (-z+zInside)**2 ) )

def Elementary_field_function_z_Mirror(j_layer, z, zInside, ro_of_layers_Mirror):
 return Density_of_charge_Mirror(j_layer, zInside, ro_of_layers_Mirror) *(zsignum(z,zInside) + (z-zInside)/sqrt( (R_of_layers[j_layer])**2 + (-z+zInside)**2 ) )

# строим график плотности заряда
# ситуация на оси симметрии облака
# заряд однороден по горизонтали, размыт по Гауссу вдоль z

names_of_layers=[ 'lower', 'middle', 'upper']       
list_of_colours=[ 'blue', 'orange' ,'green'] 
#list_of_styles=['--' , '.-', '-']  
z_range=range(10,13000)   

fig = plt.figure(figsize=(10,7))    
plt.rc('xtick', labelsize=18) 
plt.rc('ytick', labelsize=18)     
for j_layer in range (0, len(names_of_layers)):
    y=[Density_of_charge(j_layer,z, ro_of_layers)  for z in z_range];
    plt.plot(list(z_range),y, '-', label=names_of_layers[j_layer], color=list_of_colours[j_layer], linewidth=3)
#    y=[Density_of_charge_Mirror(j_layer,z, ro_of_layers_Mirror)  for z in z_range];
#    plt.plot(list(z_range),y, '--', label=names_of_layers[j_layer], color=list_of_colours[j_layer], linewidth=3)
plt.xlabel(r'altitude above the ground,$m$', rotation='horizontal', fontsize=20)
plt.ylabel(r'charge density,${C}{m}^{-3}$', rotation='vertical', fontsize=20, horizontalalignment='center', verticalalignment='bottom')
plt.legend(fontsize=20,loc=3 )
plt.show()

 # определим функцию, возвращающую поле по данным j_layer, z(высота)
     # с использованием вышезаданной плотности заряда
#      electric field sstrength in kV/m
def Field_in_z_function(z, ro_of_layers, ro_of_layers_Mirror):
    field_z=0;    
    for j_layer in range (0, len(names_of_layers)):
        def local_field_z(zInside):
            return Elementary_field_function_z(j_layer, z, zInside, ro_of_layers)    
        
        def local_field_z_Mirror(zInside):
            return Elementary_field_function_z_Mirror(j_layer, z, zInside, ro_of_layers_Mirror)    
           
        integral1, err = quad( local_field_z, zCentre_of_layers[j_layer]-zLimit,zCentre_of_layers[j_layer]+zLimit )
        integral2, err = quad( local_field_z_Mirror, zCentre_of_layers_Mirror[j_layer]-zLimit,zCentre_of_layers[j_layer]+zLimit )
        field_z+=(integral1 + integral2);
    return (10**(-3)*field_z/2/epsilon0)



field_critical_negative=[];
for j_z in range(0,  len(z_range)):
    z = z_range[j_z]
    field_critical_negative.append(10**(-3)*(-2.76e5)*0.87**((z+MountainHeight)/1000))
    


field_profile = []
for current_z in z_range:
    field_profile.append(Field_in_z_function(current_z, ro_of_layers, ro_of_layers_Mirror))
    

fig = plt.figure(figsize=(18,10)) 

# =============================================================================
# z_small_range = np.arange(zmin, z_range[j_z-1]+zstep, zstep)    
# field_critical_negative_small=[];
# 
# for j_z in range(0,  len(z_small_range)):
#     field_critical_negative_small.append(field_critical_negative[j_z])    
# =============================================================================
    
plt.plot(z_range, field_profile, linewidth=3, label='field_profile')
#plt.plot(z_range, field_critical_negative, linewidth=3, label='-critical')
#plt.plot(z_range, field_critical_positive, linewidth=3, label='+critical')
plt.xlabel(r'altitude above the ground,$m$', rotation='horizontal', fontsize=20)
plt.ylabel(r'E,$\frac{kV}{m}$', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.legend(fontsize=20,loc=1)

plt.show()    
    