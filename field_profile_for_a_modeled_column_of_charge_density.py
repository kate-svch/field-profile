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
import datetime

import matplotlib as mpl
  # инициализация констант
coef_k = 9*10**(9)    
MountainHeight=3200; 
R = 1 ;  # эффективный радиус заряженных слоёв в км


dx = 1000;  #латеральный размер ячейки в м

# fraction_number = 0, 1 - correcponds to GRAUP and SNOW - correspondingly


the_time_moment = datetime.datetime(2016, 6, 11, 11, 10) 
# LET'S LOAD ALL THE ARRAYS - from files made by "model2.py"
z_array = np.load('/home/kate-svch/Thunder/Aragats_measurements/py-codes/z_and_dens_arrays/z_array_' + datetime.datetime.strftime(the_time_moment, '%Y-%m-%d_%H:00:00') +'.npy')
dz = z_array[1] - z_array[0]  # в м


name_array = ['QGRAUP', 'QSNOW']
hydrometeors_type_quantity = len(name_array)
charge_coef_array = [1,1]

density_array = [0,0]
for jjj in range(0, len(name_array)):
    name = name_array[jjj]
    density_array[jjj] = np.load('/home/kate-svch/Thunder/Aragats_measurements/py-codes/z_and_dens_arrays/dens_array_' + name + '_'+ datetime.datetime.strftime(the_time_moment, '%Y-%m-%d_%H:00:00') +'.npy')    


def zsignum(z, zInside):
    return 2*((-z+zInside)>0)-1;

# расстояние подставляется в километрах, результата получается в kV//m (отсюда множитель 10**(-9))
# умножение на объём ячейки - чобы из плотности заряда получить собственно заряд, измеряемый в C    
def Elementary_field_function_z(fraction_number, z_index, z_Inside_index):
    return 10**(-9)*dx*dx*dz*coef_k*charge_coef_array[fraction_number]*density_array[fraction_number][z_Inside_index] *(zsignum(z_array[z_index],z_array[z_Inside_index]) + (z_array[z_index]-z_array[z_Inside_index])/sqrt( R**2 + (z_array[z_index]+z_array[z_Inside_index])**2 ) )


def Elementary_field_function_z_Mirror(fraction_number, z_index, z_Inside_index):
    return -10**(-9)*dx*dx*dz*coef_k*charge_coef_array[fraction_number]*density_array[fraction_number][z_Inside_index] *(zsignum(z_array[z_index],z_array[z_Inside_index]) + (z_array[z_index]-z_array[z_Inside_index])/sqrt( R**2 + (z_array[z_index]+z_array[z_Inside_index])**2 ) )

# строим график плотности заряда
for hydrometeor_type_ind in range (0, hydrometeors_type_quantity): 
    plt.figure(figsize=(12,4))
    plt.plot( density_array[hydrometeor_type_ind], z_array, linewidth = 3)
    plt.title('Density profile of ' + name_array[hydrometeor_type_ind] + ''+ str(the_time_moment), fontsize=22)
    plt.xlabel('density, some units', fontsize=20, horizontalalignment='right' )
    plt.ylabel('z, km', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
    plt.axis('normal')    
    plt.show()    
    
    
plt.figure(figsize=(12,4))
plt.title('Charge-density profile ' + str(the_time_moment), fontsize=22)
for hydrometeor_type_ind in range (0, hydrometeors_type_quantity): 
    plt.plot( charge_coef_array[hydrometeor_type_ind]*density_array[hydrometeor_type_ind], z_array, linewidth = 3, label = name_array[hydrometeor_type_ind])   
plt.xlabel('charge-density, '+ r'$\frac{c}{m^3}$', fontsize=20, horizontalalignment='right' )
plt.ylabel('z, km', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.axis('normal')
plt.legend(fontsize=20,loc=1)    
plt.show()      


# let's look at the "elementary_field_function"'s work:
elem_field_func_z_result = [];  # z_array_one_part = []
z_Inside_index = 150;  
fraction_number = 0;
for z_index in range(0, len(z_array)):
#for z_index in range(0, z_Inside_index):
#    z_array_one_part.append(z_array[z_index])
    if (z_Inside_index != z_index):
        elem_field_func_z_result.append(Elementary_field_function_z(fraction_number, z_index, z_Inside_index))
    else:    
        elem_field_func_z_result.append(0)
    

plt.figure(figsize=(12,4))
plt.title('Elementary_field_function_z(' + str(name_array[fraction_number]) + ', z_index, ' + str(z_Inside_index) + ' ' + str(the_time_moment), fontsize=22)
plt.plot(elem_field_func_z_result,  z_array, linewidth = 3, label = name_array[fraction_number])   
plt.xlabel('electric_field, '+ r'$\frac{kV}{m}$', fontsize=20, horizontalalignment='right' )
plt.ylabel('z, km', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.axis('normal')
plt.legend(fontsize=20,loc=1)    
plt.show()  
    

plt.figure(figsize=(12,4))
plt.title('Elementary_field_function_z(' + str(name_array[fraction_number]) + ', z_index, ' + str(z_Inside_index) + ' ' + str(the_time_moment), fontsize=22)
plt.plot(z_array, elem_field_func_z_result, linewidth = 3, label = name_array[fraction_number])   
plt.ylabel('electric_field, '+ r'$\frac{kV}{m}$', fontsize=20, horizontalalignment='right' )
plt.xlabel('z, km', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.axis('normal')
plt.legend(fontsize=20,loc=1)    
plt.show()  
    



 # определим функцию, возвращающую поле по данным на данной высоте (заданной индексом) - от всех типов гидрометеоров, с учётом отражения
     # с использованием вышезаданной плотности заряда
#      electric field sstrength in kV/m
def Field_in_z_function(z_index, density_array):
    field_z=0;    
    for bbb in range (0, hydrometeors_type_quantity ):
        for z_Inside_index in range (0, len(z_array)): 
            if (z_index != z_Inside_index):
                field_z +=  Elementary_field_function_z(bbb, z_index, z_Inside_index)
                field_z +=  Elementary_field_function_z_Mirror(bbb,  z_index, z_Inside_index)
    return field_z;                


field_critical_negative=[];
field_profile = [];
for z_ind in range(0, len(z_array)):
    current_z = z_array[z_ind]
    field_critical_negative.append(10**(-3)*(-2.76e5)*0.87**((current_z + MountainHeight)/1000))
    field_profile.append(Field_in_z_function(z_ind, density_array))
    
   

fig = plt.figure(figsize=(18,10)) 

    
plt.plot(z_array, field_profile, linewidth=3, label='field_profile')
plt.plot(z_array, field_critical_negative, linewidth=3, label='-critical')
#plt.plot(z_range, field_critical_positive, linewidth=3, label='+critical')
plt.ylabel(r'$\frac{kV}{m}$', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.legend(fontsize=20,loc=1)
plt.show()    
    