#%%  # научимся писать "предохранитель",
# позволяющий не считать и не отображать потоки больше некоторого

start=time.clock();  

ro_l0_min=1; ro_l0_max=100; 
ro_l0_quant=100;
ro_l0_range = np.linspace(ro_l0_min, ro_l0_max, ro_l0_quant) 

ro_m0_min=1; ro_m0_max=100;
ro_m0_quant=100;
ro_m0_range = np.linspace(ro_m0_min, ro_m0_max, ro_m0_quant)

flux_qq=[]; 

# чтобы слишком большие потоки не рассчитывались:
whether_we_count_the_flux=1;

#if (1==1):  
for j_ro_m0 in range(0,  len(ro_m0_range)):
    ro_m0 = ro_m0_range[j_ro_m0];
    
    flux_qq.append([]);  
      
    whether_we_count_the_flux=1;  flux=0;
    
    for j_ro_l0 in range(0,  len(ro_l0_range)):
        ro_l0 = ro_l0_range[j_ro_l0]; 
        
#        flux=0;       ЗЛО В ИСХ. КОДЕ
                   
 # предположим, что лавина только что кончилась
# если так, то у нас есть всё для нахождения потока      
        if (1<2):           
          
  # если последнее значение flux достаточно велико, предохранитель сработает сейчас              
            if (flux>12000): 
                whether_we_count_the_flux=0;
 
            if (whether_we_count_the_flux==1): 
                flux=ro_m0**2 + ro_l0**2;   
 #               flux=ro_m0*10 + ro_l0; 
                flux_qq[j_ro_m0].append(flux);
            else:
                flux_qq[j_ro_m0].append(0);

        print ("flux = " + str(flux) + ' 1/s')
    
fig = plt.figure(figsize=(12,6)) 
picture2=plt.contourf(ro_l0_range, ro_m0_range, flux_qq, 12, cmap='jet')
plt.colorbar(picture2) 
plt.title('Particles flux', fontsize=22)
plt.xlabel(r'${\rho}_{l0}$', fontsize=20, horizontalalignment='right' )
plt.ylabel(r'${\rho}_{m0}$', rotation='horizontal', fontsize=20, horizontalalignment='right', verticalalignment='top')
plt.axis('image')


elapsed=time.clock() - start
print(elapsed) 