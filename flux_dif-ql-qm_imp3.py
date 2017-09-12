from math import*
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import time

# для построения "3D"-графиков:
#%matplotlib qt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

#%%   # инициализация констант
epsilon0 = 8.85*10**(-12);  zLimit = 1000;
H=500; hm=6500-4200; hu = 9800-4200;
ql=11; qm=-60; qu=40;
sigma_l=700; sigma_m=550; sigma_u=600;
Rl=1800; Rm=3000; Ru=4000;

ro_l0=2.6*10**(-9); ro_m0=-3.3*10**(-9); ro_u0=0.7*10**(-9); 
ro_of_layers=[ ro_l0, ro_m0, ro_u0] 
ro_of_layers_Mirror=[ -ro_l0, -ro_m0, -ro_u0];    
sigma_of_layers=[ sigma_l, sigma_m, sigma_u]   
zCentre_of_layers=[H, H+hm, H+hu] 
zCentre_of_layers_Mirror=[-H, -(H+hm), -(H+hu)];
R_of_layers=[Rl, Rm, Ru]

#%%   # задаём распределение плотности зарядов

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

#%%  # строим график плотности заряда

# ситуация на оси симметрии облака
# заряд однороден по горизонтали, размыт по Гауссу вдоль z

names_of_layers=[ 'lower', 'middle', 'upper']       
list_of_colours=[ 'blue', 'orange' ,'green'] 
#list_of_styles=['--' , '.-', '-']  
zrange=range(-11000,11000)   

fig = plt.figure(figsize=(10,7))    
    
for j_layer in range (0, len(names_of_layers)):
    y=[Density_of_charge(j_layer,z, ro_of_layers)  for z in zrange];
    plt.plot(list(zrange),y, '-', label=names_of_layers[j_layer], color=list_of_colours[j_layer], linewidth=3)
    y=[Density_of_charge_Mirror(j_layer,z, ro_of_layers_Mirror)  for z in zrange];
    plt.plot(list(zrange),y, '--', label=names_of_layers[j_layer], color=list_of_colours[j_layer], linewidth=3)
plt.legend(fontsize=20,loc=1 )
plt.show()

#%%  # определим функцию, возвращающую поле по данным j_layer, z(высота)
     # с использованием вышезаданной плотности заряда
     
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
    return (field_z/2/epsilon0)

#%%  # графики зависимости потока и длины лавины от зарядов нижнего и среднего слоёв,
    # " 2D с заливкой "
    # добавим график поля на поверхности, построенный в тех же осях

#   ОПТИМИЗИРОВАНО:
#   когда лавина кончилась, можно дальше профиль поля не считать
#   критическое поле negative задаётся сразу и навсегда, positive не требуется
#   написать поиск начала и конца лавины дихотомией
#   нахождение поля на данной высоте вынесено в отдельную функцию

#   ОПТИМИЗИРОВАТЬ?:   
#   включить память предыдущих начала и конца ?
 
# ситуация на оси симметрии облака
# заряд однороден по горизонтали, размыт по Гауссу вдоль z

zmin=0; zmax=5000; zstep=10;   MountainHeight=1000;  
# при резкой смене значений зарядов или высот убедиться,
# что нет отрицательной лавины за пределами (zmin, zmax)
quantityofiterations = 8;
z_range = np.arange(zmin, zmax + zstep, zstep) 
names_of_layers=[ 'lower', 'middle', 'upper']       
list_of_colours=[ 'blue', 'orange' ,'green'] 
list_of_styles=['--' , '.-', '-']  

start=time.clock();  

ro_l0_min=2.2e-9; ro_l0_max=2.55e-9; 
ro_l0_quant=16;
ro_l0_range = np.linspace(ro_l0_min, ro_l0_max, ro_l0_quant) 

ro_m0_min=-3.2e-9; ro_m0_max=-3.55e-9;
ro_m0_quant=16;
ro_m0_range = np.linspace(ro_m0_min, ro_m0_max, ro_m0_quant)

field_critical_negative=[];  avlength_qq=[];  flux_qq=[]; field_ground_qq=[]; 
field_deriv_qq=[];  # будем писать ещё и производную поля в конце лавины


field_critical_negative=[];
for j_z in range(0,  len(z_range)):
    z = z_range[j_z]
    field_critical_negative.append(-2.76e5*0.87**((z+MountainHeight)/1000))

#if (1==1):  
for j_ro_m0 in range(0,  len(ro_m0_range)):
    ro_m0 = ro_m0_range[j_ro_m0];
    ro_of_layers[1]=ro_m0;
    ro_of_layers_Mirror[1]=-ro_m0;
    
    avlength_qq.append([]);  flux_qq.append([]);   field_ground_qq.append([]); field_deriv_qq.append([]);
    
    # чтобы слишком большие потоки не рассчитывались:
    whether_we_count_the_flux=1;    flux=1;
    
    for j_ro_l0 in range(0,  len(ro_l0_range)):
        ro_l0 = ro_l0_range[j_ro_l0]; 
        ro_of_layers[0]=ro_l0;
        ro_of_layers_Mirror[0]=-ro_l0;
        
        avalanche_length=0; avalanche_indicator=0; field_derivative=0;
        field_profile=[];         
            
 #  это - нахождение потока без "предохранителя ужасной бомбардировки частицами"      
 #       flux=exp( avalanche_length**2/2/7.3e6 * field_derivative ) 
        
        if (flux>1.0e2): 
            whether_we_count_the_flux=0;
  
        
        if (whether_we_count_the_flux==1): 
            
            j_z=0;    z = z_range[j_z];   avalanche_end =avalanche_start =0;
            
            field_ground_qq[j_ro_m0].append(Field_in_z_function(z, ro_of_layers, ro_of_layers_Mirror))
            
            while (avalanche_indicator<2)and(j_z<len(z_range)):
                z = z_range[j_z]
                field_profile.append(Field_in_z_function(z, ro_of_layers, ro_of_layers_Mirror))
           
        # нашли поле в данной точке. Посмотрим, а не в лавине ли мы:
                 
                if (field_profile[j_z] < field_critical_negative[j_z])and (avalanche_indicator<1):
        #            avalanche_length+=1;   # длина лавины в единицах zstep
                    avalanche_indicator=1;   # равен нулю, если лавины ещё не было
                    
                    rightedge=z_range[j_z];   leftedge=z_range[j_z-1];
                    avalanche_start = (leftedge + rightedge)/2
     #               quantityofiterations = 8;
                
                    for q_to_find_the_precise_av_start in range(0, quantityofiterations):
                        
                        z = avalanche_start;  
                        field_current=(Field_in_z_function(z, ro_of_layers, ro_of_layers_Mirror));
                        field_crit_negat_current = -2.76e5*0.87**((z+MountainHeight)/1000);
                        
                        if field_current > field_crit_negat_current : leftedge = avalanche_start
                        else:   rightedge = avalanche_start
                        avalanche_start = (leftedge + rightedge)/2
    
                 # мы сейчас в нижней точке лавины: найдём производную в её конце           
                    
                    field_left=Field_in_z_function(leftedge, ro_of_layers,ro_of_layers_Mirror);
                    field_right=Field_in_z_function(rightedge, ro_of_layers,ro_of_layers_Mirror);
            
                    field_derivative=abs((field_right-field_left))/(rightedge-leftedge); 
                        
             # проверим, - а вдруг лавина только что кончилась?      
                if (field_profile[j_z] > field_critical_negative[j_z])and (avalanche_indicator==1):           
                    avalanche_indicator=2;  # 1 в лавине, 2 - если уже прошла и пора сворачиваться
                    rightedge=z_range[j_z];   leftedge=z_range[j_z-1];
                    avalanche_end = (leftedge + rightedge)/2
    #                quantityofiterations = 8;
                
                    for q_to_find_the_precise_av_start in range(0, quantityofiterations):
                        
                        z = avalanche_end;  
                        field_current=(Field_in_z_function(z, ro_of_layers, ro_of_layers_Mirror));
                        field_crit_negat_current = -2.76e5*0.87**((z+MountainHeight)/1000)
                                                      
                        if field_current < field_crit_negat_current : leftedge = avalanche_end
                        else:   rightedge = avalanche_end                
                        avalanche_end = (leftedge + rightedge)/2
                        
                    avalanche_length = avalanche_end - avalanche_start; 
                    avlength_qq[j_ro_m0].append(avalanche_length);
                
    #                flux=exp( avalanche_length**2/2/7.3e6 * field_derivative ) 
     #               print ("I've just counted the flux value" )
            
                j_z+=1;
# будем строить каждый третитй профиль поля, чтобы убедиться в правдоподобности итоговых результатов                       
            if ((j_ro_m0+1)%4==0)&((j_ro_l0+1)%4==0):
                fig = plt.figure(figsize=(16,8)) 
                z_small_range = np.arange(zmin, z_range[j_z-1]+zstep, zstep)    
                field_critical_negative_small=[];
                
                for j_z in range(0,  len(z_small_range)):
                    field_critical_negative_small.append(field_critical_negative[j_z])    
                plt.title(r'${\rho}_{m0} = $'+str(ro_m0)+r'$; {\rho}_{l0} = $'+str(ro_l0)+r'$; field-der = $'+str(field_derivative), fontsize=16)    
                plt.plot(z_small_range, field_profile, linewidth=3, label='field_profile')
                plt.plot(z_small_range, field_critical_negative_small, linewidth=3, label='-critical')
                plt.legend(fontsize=20,loc=1)
                plt.show()
            
            flux=exp( avalanche_length**2/2/7.3e6 * field_derivative )    
            flux_qq[j_ro_m0].append(flux);
            field_deriv_qq[j_ro_m0].append(field_derivative);
        else:
            flux_qq[j_ro_m0].append(1);
            avlength_qq[j_ro_m0].append(1);
            field_ground_qq[j_ro_m0].append(1);
            field_deriv_qq[j_ro_m0].append(1);
 
#        print ("I've just counted the flux value" )    
#        print ("length of avalanche = " + str(avalanche_length) + ' m')
#        print ("flux = " + str(flux) + ' 1/s')
    
fig = plt.figure(figsize=(8,6)) 
picture1=plt.contourf(ro_l0_range, ro_m0_range, avlength_qq, 12, cmap='jet')
plt.colorbar(picture1) 
plt.title('Avalanche length', fontsize=22)
plt.xlabel(r'${\rho}_{l0}$', fontsize=20, horizontalalignment='right' )
plt.ylabel(r'${\rho}_{m0}$', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.axis('image')

fig = plt.figure(figsize=(8,6)) 
picture2=plt.contourf(ro_l0_range, ro_m0_range, field_deriv_qq, 12, cmap='jet')
plt.colorbar(picture2) 
plt.title('Field derivative', fontsize=22)
plt.xlabel(r'${\rho}_{l0}$', fontsize=20, horizontalalignment='right' )
plt.ylabel(r'${\rho}_{m0}$', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.axis('image')

fig = plt.figure(figsize=(8,6)) 
picture2=plt.contourf(ro_l0_range, ro_m0_range, flux_qq, 12, cmap='jet')
plt.colorbar(picture2) 
plt.title('Particles flux', fontsize=22)
plt.xlabel(r'${\rho}_{l0}$', fontsize=20, horizontalalignment='right' )
plt.ylabel(r'${\rho}_{m0}$', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.axis('image')
    
fig = plt.figure(figsize=(8,6)) 
picture3=plt.contourf(ro_l0_range, ro_m0_range, field_ground_qq, 12, cmap='jet')
plt.colorbar(picture3) 
plt.title('Field strength near the ground', fontsize=22)
plt.xlabel(r'${\rho}_{l0}$', fontsize=20, horizontalalignment='right' )
plt.ylabel(r'${\rho}_{m0}$', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.axis('image')

elapsed=time.clock() - start
print(elapsed)  
                    
#%% построим зависимость производной поля от ro_l0 при максимальном ro_m0 из диапазона

plt.title(r'${\rho}_{m0} = $'+str(ro_m0) + r'$; field_der ({\rho}_{l0}) $', fontsize=16)    
plt.plot(ro_l0_range, field_deriv_qq[-1], linewidth=3)
plt.xlabel(r'${\rho}_{m0}$', fontsize=20, horizontalalignment='right' )
plt.ylabel(r'$field deriv$', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.show()



#%% попробуем заливать каждый интервал своим цветом БЕЗ РАЗМЫТИЙ НА ГРАНИЦАХ
# для длины лавины и производной поля - границы диапазонов в levels не подобраны !!!
'''
fig = plt.figure(figsize=(8,6)) 
levels = [1, 1.1, 1.2, 1.3, 1.7]
picture1=plt.contourf(ro_l0_range, ro_m0_range, avlength_qq, levels, colors=('r', 'y','g', 'b'))
picture1.cmap.set_under('orange')
picture1.cmap.set_over('cyan')
plt.colorbar(picture1) 
plt.title('Avalanche length', fontsize=22)
plt.xlabel(r'${\rho}_{l0}$', fontsize=20, horizontalalignment='right' )
plt.ylabel(r'${\rho}_{m0}$', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.axis('image')

fig = plt.figure(figsize=(8,6)) 
levels = [1, 1.1, 1.2, 1.3, 1.7]
picture1=plt.contourf(ro_l0_range, ro_m0_range, field_deriv_qq, levels, colors=('r', 'y','g', 'b'))
picture1.cmap.set_under('orange')
picture1.cmap.set_over('cyan')
plt.colorbar(picture1) 
plt.title('Field derivative', fontsize=22)
plt.xlabel(r'${\rho}_{l0}$', fontsize=20, horizontalalignment='right' )
plt.ylabel(r'${\rho}_{m0}$', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.axis('image')   '''

fig = plt.figure(figsize=(8,6)) 
levels = [1, 1.1, 1.2, 1.3, 1.7]
picture1=plt.contourf(ro_l0_range, ro_m0_range, flux_qq, levels, colors=('r', 'y','g', 'b'))
picture1.cmap.set_under('orange')
picture1.cmap.set_over('cyan')
plt.colorbar(picture1) 
plt.title('Particles flux', fontsize=22)
plt.xlabel(r'${\rho}_{l0}$', fontsize=20, horizontalalignment='right' )
plt.ylabel(r'${\rho}_{m0}$', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.axis('image')



# of colors.

# Our data range extends outside the range of levels; make
# data below the lowest contour level yellow, and above the
# highest level cyan:


#%% вывод посчитанных выше картинок с заливкой только интересующего диапазона
# ОПРЕДЕЛИТЬ И ПОДСТАВИТЬ MAX, MIN VALUE - ПО ПРЕДЫДУЩЕМУ

quantity_of_parts=10; minimum_value=0;  maximum_value=1000;
#size_of_part=maximum_value/quantity_of_parts;
value_bounds = np.linspace(minimum_value, maximum_value, quantity_of_parts+1);

cpool=[];
for i in range(quantity_of_parts):
    cpool.append('white');
cpool[-2]='green';   
#cpool = ['red','orange','yellow','green','cyan', 'blue', 'violet', 'grey', 'black', 'red' ]
 
fig = plt.figure(figsize=(30,4)) 
#fig = plt.figure(1,(14,4))
ax = fig.add_subplot(131)    # ШАГ РЕГУЛЯРНЫЙ!! (в следующей строке)
# зачем 'indexed', неясно
cmap = mpl.colors.ListedColormap(cpool, 'indexed') # Задаём дискретную шкалу цветов из списка cpool
cs = ax.pcolor(ro_l0_range, ro_m0_range, avlength_qq, cmap=cmap) # вызываем метод pcolor. Вводим пользовательскую раскраску через cmap
#cbar = fig.colorbar(cs, ticks=value_bounds) # рисуем шкалу colorbar для изображения cs. В качестве черточек шкалы указываем bounds. 
cbar = fig.colorbar(cs)

#%%  for flux

quantity_of_parts=16; minimum_value=1.002;  maximum_value=1.02;
#size_of_part=maximum_value/quantity_of_parts;
value_bounds = np.linspace(minimum_value, maximum_value, quantity_of_parts+1);

cpool=[];
for i in range(quantity_of_parts):
    cpool.append('white');
cpool[-3]='green';   

fig = plt.figure(figsize=(30,4)) 
#fig = plt.figure(1,(14,4))
ax = fig.add_subplot(131)    # ШАГ РЕГУЛЯРНЫЙ!! (в следующей строке)
cmap = mpl.colors.ListedColormap(cpool, 'indexed') # Задаём дискретную шкалу цветов из списка cpool
cs = ax.pcolor(ro_l0_range, ro_m0_range, flux_qq, cmap=cmap) # вызываем метод pcolor. Вводим пользовательскую раскраску через cmap
#cbar = fig.colorbar(cs, ticks=value_bounds) # рисуем шкалу colorbar для изображения cs. В качестве черточек шкалы указываем bounds. 
cbar = fig.colorbar(cs)

#%%  for field

quantity_of_parts=10; minimum_value=30000;  maximum_value=55000;
#size_of_part=maximum_value/quantity_of_parts;
value_bounds = np.linspace(minimum_value, maximum_value, quantity_of_parts+1);

cpool=[];
for i in range(quantity_of_parts):
    cpool.append('white');
cpool[-3]='green';   

fig = plt.figure(figsize=(30,4)) 
#fig = plt.figure(1,(14,4))
ax = fig.add_subplot(131)    # ШАГ РЕГУЛЯРНЫЙ!! (в следующей строке)
cmap = mpl.colors.ListedColormap(cpool, 'indexed') # Задаём дискретную шкалу цветов из списка cpool
cs = ax.pcolor(ro_l0_range, ro_m0_range, field_ground_qq, cmap=cmap) # вызываем метод pcolor. Вводим пользовательскую раскраску через cmap
#cbar = fig.colorbar(cs, ticks=value_bounds) # рисуем шкалу colorbar для изображения cs. В качестве черточек шкалы указываем bounds. 
cbar = fig.colorbar(cs)

#%%  тест автоматического задания раскраски "only one"
cpool=[];
for i in range(quantity_of_parts):
    cpool.append('white');

cpool[-3]='green';    

#%% вывод последних графиков длины лавины и потока

fig = plt.figure(figsize=(8,6))     
plt.plot(ro_l0_range, avlength_qq[j_ro_m0], linewidth=3, label='av-length', color='orange')
plt.show()

fig_flux = plt.figure(figsize=(8,6))     
plt.plot(ro_l0_range, flux_qq[j_ro_m0], linewidth=3, label='av-length', color='violet')
plt.show()

#%%  # строим профиль поля для последней плотности, ВСЕ СЛОИ ВМЕСТЕ С ОТРАЖЕНИЯМИ 
fig = plt.figure(figsize=(8,6)) 

z_small_range = np.arange(zmin, z_range[j_z-1]+zstep, zstep)    
field_critical_negative_small=[];

for j_z in range(0,  len(z_small_range)):
    field_critical_negative_small.append(field_critical_negative[j_z])    
    
plt.plot(z_small_range, field_profile, linewidth=3, label='field_profile')
plt.plot(z_small_range, field_critical_negative_small, linewidth=3, label='-critical')
#plt.plot(z_range, field_critical_positive, linewidth=3, label='+critical')
plt.legend(fontsize=20,loc=1)

plt.show()

