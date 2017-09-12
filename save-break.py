#%%  # научимся писать "предохранитель",
# позволяющий не считать и не отображать потоки больше некоторого

from math import*
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import time

# для построения "3D"-графиков:
#%matplotlib qt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

start=time.clock();  

x_min=1; x_max=100; 
x_quant=100;
x_range = np.linspace(x_min, x_max, x_quant) 

y_min=1; y_max=100;
y_quant=100;
y_range = np.linspace(y_min, y_max, y_quant)

test_qq=[]; 

# чтобы слишком большие потоки не рассчитывались:
whether_we_count_the_flux=1;

#if (1==1):  
for j_y in range(0,  len(y_range)):
    y = y_range[j_y];
    
    test_qq.append([]);  
      
    whether_we_count_the_flux=1;  flux=0;
    
    for j_x in range(0,  len(x_range)):
        x = x_range[j_x]; 
        
#        flux=0;       ЗЛО В ИСХ. КОДЕ
                   
 # предположим, что лавина только что кончилась
# если так, то у нас есть всё для нахождения потока      
        if (1<2):           
          
  # если последнее значение flux достаточно велико, предохранитель сработает сейчас              
            if (flux>12000): 
                whether_we_count_the_flux=0;
 
            if (whether_we_count_the_flux==1): 
                flux=y**2 + x**2;   
 #               flux=ro_m0*10 + ro_l0; 
                test_qq[j_y].append(flux);
            else:
                test_qq[j_y].append(0);

#        print ("flux = " + str(flux) + ' 1/s')
    
fig = plt.figure(figsize=(12,6)) 
picture2=plt.contourf(x_range, y_range, test_qq, 12, cmap='jet')
plt.colorbar(picture2) 
plt.title('Test', fontsize=22)
plt.xlabel(r'${\rho}_{l0}$', fontsize=20, horizontalalignment='right' )
plt.ylabel(r'${\rho}_{m0}$', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.axis('image')


elapsed=time.clock() - start
print(elapsed) 

#%%

quantity_of_parts=10; minimum_value=0;  maximum_value=13000;
#size_of_part=maximum_value/quantity_of_parts;
value_bounds = np.linspace(minimum_value, maximum_value, quantity_of_parts+1);

cpool=[];
for i in range(quantity_of_parts):
    cpool.append('white');
cpool[-3]='green';   

fig = plt.figure(figsize=(20,5)) 
#fig = plt.figure(1,(14,4))
ax = fig.add_subplot(131)    # ШАГ РЕГУЛЯРНЫЙ!! (в следующей строке)
cmap = mpl.colors.ListedColormap(cpool, 'indexed') # Задаём дискретную шкалу цветов из списка cpool
cs = ax.pcolor(x_range, y_range, test_qq, cmap=cmap) # вызываем метод pcolor. Вводим пользовательскую раскраску через cmap
#cbar = fig.colorbar(cs, ticks=value_bounds) # рисуем шкалу colorbar для изображения cs. В качестве черточек шкалы указываем bounds. 
cbar = fig.colorbar(cs)

