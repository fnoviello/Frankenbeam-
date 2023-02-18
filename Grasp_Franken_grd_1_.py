#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 20:22:58 2022

@author: fabionoviello
"""

'''

          READ GRASP phi theta, extract cuts ad a perturbation
          to a specifc theta angular scale and repack into GRASP .grd format.
          
          
 
                     Fabio Noviello
 
                         (2022-2023)
 
 ''' 
 
 ### Libraries
 
 
 
import numpy as np

import matplotlib.pyplot as plt



#--parameters

h_rowskip = 12  # Header lines to skip

npoints_theta = 901

ncuts_phi = 36 # phi grid points (see header file) # the "cuts" are a legacy from cuts-oriented code.

cut_space = 360/ncuts_phi # spacing between phi grid points (if evenly spaced)

#The grasp grid is ordered as phi-theta, with phi growing faster than theta,
# e.g. (phi, theta) = (0, 0), (1,0), (2,0).... However, I want to extract a 
# (theta, phi) cut such as (0,0)), (1, 0)..., so in grasp (phi theta) notation this reads
# as (0,0), (0,1), (0,2)...etc. To achieve this I must read (phi = 0, theta = 0) in the grasp file,
# skip n_cuts_phi rows and read the ncuts_theta + 1 row, which gives me (phi = 0, theta = 1)
# and so on.. therefor:


phi_rowskip = ncuts_phi 



ntheta = 901 # theta points in each GRASP phi cut (S)
#ntheta_half = int(ntheta/2 + 1) #halfway/peak point SYMMETRICAL BEAMS ONLY.

theta_min = 0.0          # Grid extrema (file dependent)
theta_max = 90.0                  ##
phi_min = 0.0                     ##
phi_max = 360.0                   ##


x_tick_min = 0 # min for x_tick axis
x_tick_max = 90 # min for x_tick axis

########################################################
#--more parameters-this time related to dB increment in cuts--
######################################################


#--test with same increment for all cuts-nomenclature already set for 1st cut 
#out of X -will need to expamd.. 


#--THE FOLLOWING PARAMETERS WILL BE COMMON FOR ALL CUTS:
    
theta_range = 90 #theta range (in degrees) for each cut (see beam header file)

theta_fac = ntheta/theta_range # angles to data points conversion factor    
                               # i.e. how many data points for each angle ?                     



delta_db_1 =  60# By how much we want to increase) sidelobes (in dB)

strd_db_1 = str(delta_db_1) #for plots etc..



delta_1 = 10**(delta_db_1/10) # linear version of the above, assuming that
                          # the db pattern is defined as 10*np.log10(linear)
                          
modtheta_min_1_1 = 30  # initial theta for Frankenbeam enhancement if choosing theta
                # range option (see code below).
                # 1_1 stands for "1st min, 1st cut" etc.
modtheta_max_1_1 = 60 # final theta for Frankenbeam enhancement if choosing theta
                # range option    


modtheta_min_1_1 = int(modtheta_min_1_1*theta_fac) #min and max of modified
modtheta_max_1_1 = int(modtheta_max_1_1*theta_fac) # interval counted as 
                                                     #data points



#---------------
# ------Read data

infile = "/Users/fabionoviello/Documents/Research/Litebird//Beams_Cristian_2022_03/BB1_theta_phi/los.grd"


indata = np.loadtxt(infile, skiprows = h_rowskip, unpack = True)



print("Filename: ", infile)
print('          ')


insiz = np.shape(indata)
print("Input data (minus heade+r) size is  ", insiz)


#Grasp conventions 

ReCo = indata[0, :] # GRASP polarised beams conventions.
ImCo = indata[1, :]
ReCx = indata[2, :]
ImCx = indata[3, :]

CoBeam = (ReCo**2) + (ImCo**2)   # Power of the above
CxBeam = (ReCx**2) + (ImCx**2)


Cosiz = np.shape(CoBeam)
print("Copolar beam fle size is  ", Cosiz)
Cxsiz = np.shape(CxBeam)
print("Crosspolar beam fle size is  ", Cxsiz)

CoBeamgrd = CoBeam.reshape((ncuts_phi, npoints_theta))
CxBeamgrd = CxBeam.reshape((ncuts_phi, npoints_theta))
CoBeamgrd_dB=10*np.log10(CoBeamgrd)
CxBeamgrd_dB=10*np.log10(CxBeamgrd)

print("max copolar beam is ", np.max(CoBeamgrd))
print("min copolar beam is ", np.min(CoBeamgrd))
print("max crosspolar beam is ", np.max(CxBeamgrd))
print("min crosspolar beam is ", np.min(CxBeamgrd))

print("max copolar beam (dB)is ", np.max(CoBeamgrd_dB))
print("min copolar beam (dB) is ", np.min(CoBeamgrd_dB))
print("max crosspolar beam (dB) is ", np.max(CxBeamgrd_dB))
print("min crosspolar beam (dB) is ", np.min(CxBeamgrd_dB))



##--NOW THE LOOP

# Define matriXes containing all cuts



ReCo_Mat = np.zeros((ncuts_phi,npoints_theta)) 
ImCo_Mat = np.zeros((ncuts_phi,npoints_theta)) 
PhasesCo_Mat = np.zeros((ncuts_phi,npoints_theta)) #Phase of Re wavefront @this_point

ReCx_Mat = np.zeros((ncuts_phi,npoints_theta)) 
ImCx_Mat = np.zeros((ncuts_phi,npoints_theta)) 
PhasesCx_Mat = np.zeros((ncuts_phi,npoints_theta)) #Phase of Cx wavefront @this_point

CoBeam_Mat_Pwr = np.zeros((ncuts_phi,npoints_theta)) # Power
CxBeam_Mat_Pwr = np.zeros((ncuts_phi,npoints_theta)) # Power
CoFrank_Mat_Pwr = np.zeros((ncuts_phi,npoints_theta)) # Franken (modified) beam Power
CxFrank_Mat_Pwr = np.zeros((ncuts_phi,npoints_theta)) #  "    "   
CoFrank_Mat_Pwr_dB=10*np.log10(CoFrank_Mat_Pwr)
CxFrank_Mat_Pwr_dB=10*np.log10(CxFrank_Mat_Pwr)

Re_CoFrank_Mat = np.zeros((ncuts_phi,npoints_theta)) # Real part of Co Frankenbeam 
                                                     # for Grasp.grd output
Im_CoFrank_Mat = np.zeros((ncuts_phi,npoints_theta)) # Imaginary part of  Co Frankenbeam 
                                                     # for Grasp.grd output
Re_CxFrank_Mat = np.zeros((ncuts_phi,npoints_theta)) # Real part of Cx Frankenbeam 
                                                     # for Grasp.grd output
Im_CxFrank_Mat = np.zeros((ncuts_phi,npoints_theta)) # Imaginary part of Cx Frankenbeam 
                                                     # for Grasp.grd output
                                                     
Re_CoFrank_out = np.zeros((ncuts_phi*npoints_theta))
Im_CoFrank_out = np.zeros((ncuts_phi*npoints_theta))
Re_CxFrank_out = np.zeros((ncuts_phi*npoints_theta))       
Im_CxFrank_out = np.zeros((ncuts_phi*npoints_theta)) 
        
for h in range(0, ncuts_phi):
    for i in range(0, npoints_theta):  
       
            
           # j is also incremented by h for every cut (see file descriptor.)
           
            j = i*ncuts_phi + h     #should be i*rowskip_phi to be fair, but the latter is the same as 
                                         # phi rowskip from the parameters section....
            
            ReCo_Mat[h, i] = ReCo[j]
            ImCo_Mat[h,i] = ImCo[j]
            PhasesCo_Mat[h,i] = np.arctan(ImCo[j]/ReCo[j]) # in radians
            ReCx_Mat[h, i] = ReCx[j]
            ImCx_Mat[h,i] = ImCx[j]
            PhasesCx_Mat[h,i] = np.arctan(ImCx[j]/ReCx[j]) # in radians
            
            
            
            
    PlotCo = ReCo_Mat[h,:]**2 + ImCo_Mat[h,:]**2
    PlotCx = ReCx_Mat[h,:]**2 + ImCx_Mat[h,:]**2
    
    PlotCo_dB = 10*np.log10(PlotCo)
    PlotCx_dB = 10*np.log10(PlotCx)
         
    #   print('h is  ', h, "i is ", i, "j is ", j)
          #  print("i is ", i)
          #  print("j is ", j)
       
      #  FROM HERE again NICE LOOPY PLOTS
            
    x_axis = np.linspace(0, npoints_theta, num = npoints_theta)
    x_tickaxis = np.linspace(0, npoints_theta, 10)
    x_ticklabels = (np.linspace (0, 90, 11))
    x_ticklabels = ('0', '10', '20', '30', '40', \
    '50', '60', '70', '80', '90') #change according to your beam.
# x_ticklabels = ('-100', '-60', '-20','0', 
#                 '20', '60', '100') #change according to your beam.
  
    fig =plt.figure()
   
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])            
    ax.plot(x_axis,PlotCo_dB, color = 'b', label = "Co")
    ax.plot(x_axis,PlotCx_dB, color = 'r', label = "Cx")
    ax.set_xlabel('Theta (deg)')
    ax.set_ylabel('dB')
    ax.set_title(' Input cut phi = ' + str(int(h*cut_space)))# for integer cut spacing
   # ax.set_title(' Input cut phi = ' + str(h*cut_space)) # otherwise
    ax.set_xticks(x_tickaxis)
    ax.set_xticklabels(x_ticklabels)
    plt.legend(loc="upper right")
    plt.show(block=False)

##################################################################
##TLOOP ALTERATION STARTS HERE
##
## APPLYING THE SAME DB INCREASE TO EACH CUT
##
## This can be modified to add different incresases to different cuts
## in different theta ranges and only co or cx if required.
##
## To achieve the above the initial parameters can be modified and the 
## lopp sequence fine-tuned.
##
##################################################################




for h in range(0, ncuts_phi):
    for i in range(0, npoints_theta):  
       
            
           # j is also incremented by h for every cut (see file descriptor.)
           
            j = i*ncuts_phi + h 
            
            ReCo_Mat[h, i] = ReCo[j]
            ImCo_Mat[h,i] = ImCo[j]
            ReCx_Mat[h, i] = ReCx[j]
            ImCx_Mat[h,i] = ImCx[j]

 



            CoBeam_Mat_Pwr[h,i] = (ReCo_Mat[h,i]**2) + (ImCo_Mat[h,i]**2)   # GRASP polarised beams conventions.
            CxBeam_Mat_Pwr[h,i] = (ReCx_Mat[h,i]**2) + (ImCx_Mat[h,i]**2) 
            
                        
    CoFrank_Mat_Pwr[h] = CoBeam_Mat_Pwr[h] # first using Cobeam as a template...
    CxFrank_Mat_Pwr[h] = CxBeam_Mat_Pwr[h] # first using Cxbeam as a template...



    CoFrank_Mat_Pwr[h, modtheta_min_1_1 : modtheta_max_1_1] = CoBeam_Mat_Pwr[h, modtheta_min_1_1 : modtheta_max_1_1]* delta_1
    CxFrank_Mat_Pwr[h, modtheta_min_1_1 : modtheta_max_1_1] = CxBeam_Mat_Pwr[h, modtheta_min_1_1 : modtheta_max_1_1]* delta_1


    CoFrank_Mat_Pwr_dB[h] = 10*np.log10(CoFrank_Mat_Pwr[h])
    CxFrank_Mat_Pwr_dB[h] = 10*np.log10(CxFrank_Mat_Pwr[h])
#
 
    x_axis = np.linspace(0, npoints_theta, num = npoints_theta)
    x_tickaxis = np.linspace(0, npoints_theta, 10)
    x_ticklabels = (np.linspace (0, 90, 11))
    x_ticklabels = ('0', '10', '20', '30', '40', \
    '50', '60', '70', '80', '90') #change according to your beam.
# x_ticklabels = ('-100', '-60', '-20','0', 
#                 '20', '60', '100') #change according to your beam.
  
    fig =plt.figure()
   
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])            
    ax.plot(x_axis,CoFrank_Mat_Pwr_dB[h], color = 'b', label = "Co")
    ax.plot(x_axis,CxFrank_Mat_Pwr_dB[h], color = 'r', label = "Cx")
    ax.set_xlabel('Theta (deg)')
    ax.set_ylabel('dB')
    ax.set_title(' Frankencut phi = ' + str(int(h*cut_space)))# for integer cut spacing
   # ax.set_title(' Input cut phi = ' + str(h*cut_space)) # otherwise
    ax.set_xticks(x_tickaxis)
    ax.set_xticklabels(x_ticklabels)
    plt.legend(loc="upper right")
    plt.show(block=False)     

#Creating real and imaginary parts of Co and Cx Frank for output grid
    
for h in range(0, ncuts_phi):
    for i in range(0, npoints_theta):  
    
        Re_CoFrank_Mat[h,i] = np.sqrt(CoFrank_Mat_Pwr[h,i])*np.cos(PhasesCo_Mat[h,i])    
        Im_CoFrank_Mat[h,i] = np.sqrt(CoFrank_Mat_Pwr[h,i])*np.sin(PhasesCo_Mat[h,i])                      
        Re_CxFrank_Mat[h,i] = np.sqrt(CxFrank_Mat_Pwr[h,i])*np.cos(PhasesCx_Mat[h,i])    
        Im_CxFrank_Mat[h,i] = np.sqrt(CxFrank_Mat_Pwr[h,i])*np.sin(PhasesCx_Mat[h,i]) 

        #j = i*ncuts_phi + h 
        
        # placeholder from input ReCo_Mat[h, i] = ReCo[j]
        
        j=i*ncuts_phi  + h  

        Re_CoFrank_out[j]= Re_CoFrank_Mat[h, i]  
        Im_CoFrank_out[j]= Im_CoFrank_Mat[h, i] 
        Re_CxFrank_out[j]= Re_CxFrank_Mat[h, i]  
        Im_CxFrank_out[j]= Im_CxFrank_Mat[h, i]                            
 
##############################################    
#### 
####---Now repack Mats into Grasp .grd format.
####
##############################################



with open ('/Users/fabionoviello/Desktop/TestFrankenbeam2.grd', 'a') as f :
    f.writelines(['Field data in grid', '    ', '     ',
                    '\n ++++ '
                    '\n        1           3           2           7  '
                    '\n        0           0   '
                    
                    '\n 0.0000000000E+00  0.0000000000E+00   0.3600000000E+03    0.9000000000E+02 '
                    '\n            36         901           0 \n'  ])
    #for h in range(0, ncuts_phi):
    for i in range(0, npoints_theta*ncuts_phi):  
       
          
            
           f.write("{a: 10.10e} {b: 10.10e} {c: 10.10e} {d: 10.10e}\n".format(a = Re_CoFrank_out[i],
                                                                      b= Im_CoFrank_out[i],
                                                                      c= Re_CxFrank_out[i],
                                                                      d= Im_CxFrank_out[i]))
          

f.close()   



##---------------
#
#  Now reopen the Frankenbeam Grid to check the results and display as cuts
#
#         THIS CROSSCHECK CAN BE COMMENTED OUT IF NOT REQUIRED!!
#
##---------------

h_rowskip_chk = 6  # Header lines to skip


infile_chk = "/Users/fabionoviello/Desktop/TestFrankenbeam2.grd"
   
   

indata_chk = np.loadtxt(infile_chk, skiprows = h_rowskip_chk, unpack = True)



print("Crosscheck filename: ", infile_chk)
print('          ')


insiz_chk = np.shape(indata_chk)
print("Crosscheck data (minus header) size is  ", insiz_chk)

#-- Plot the modified output grid 

ReCo_chk = indata_chk[0, :] # GRASP polarised beams conventions.
ImCo_chk = indata_chk[1, :]
ReCx_chk = indata_chk[2, :]
ImCx_chk = indata_chk[3, :]

ReCo_chk_Mat = np.zeros((ncuts_phi,npoints_theta)) 
ImCo_chk_Mat = np.zeros((ncuts_phi,npoints_theta))

ReCx_chk_Mat = np.zeros((ncuts_phi,npoints_theta)) 
ImCx_chk_Mat = np.zeros((ncuts_phi,npoints_theta)) 

for h in range(0, ncuts_phi):
    for i in range(0, npoints_theta):  
       
            
           # j is also incremented by h for every cut (see file descriptor.)
           
            j = i*ncuts_phi + h     #should be i*rowskip_phi to be fair, but the latter is the same as 
                                         # phi rowskip from the parameters section....
            
            ReCo_chk_Mat[h,i] = ReCo_chk[j]
            ImCo_chk_Mat[h,i] = ImCo_chk[j]
            ReCx_chk_Mat[h,i] = ReCx_chk[j]
            ImCx_chk_Mat[h,i] = ImCx_chk[j]
            
    PlotCo_chk = ReCo_chk_Mat[h,:]**2 + ImCo_chk_Mat[h,:]**2
    PlotCx_chk = ReCx_chk_Mat[h,:]**2 + ImCx_chk_Mat[h,:]**2
    
    PlotCo_chk_dB = 10*np.log10(PlotCo_chk)
    PlotCx_chk_dB = 10*np.log10(PlotCx_chk) 
  
    
    x_axis = np.linspace(0, npoints_theta, num = npoints_theta)
    x_tickaxis = np.linspace(0, npoints_theta, 10)
    x_ticklabels = (np.linspace (0, 90, 11))
    x_ticklabels = ('0', '10', '20', '30', '40', \
    '50', '60', '70', '80', '90') #change according to your beam.
# 
  
    fig =plt.figure()
   
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])            
    ax.plot(x_axis,PlotCo_chk_dB, color = 'b', label = "Co")
    ax.plot(x_axis,PlotCx_chk_dB, color = 'r', label = "Cx")
    ax.set_xlabel('Theta (deg)')
    ax.set_ylabel('dB')
    ax.set_title(' Frankenbeam output check phi = ' + str(int(h*cut_space)))# for integer cut spacing
   # ax.set_title(' Input cut phi = ' + str(h*cut_space)) # otherwise
    ax.set_xticks(x_tickaxis)
    ax.set_xticklabels(x_ticklabels)
    plt.legend(loc="upper right")
    plt.show(block=False)

            
            
##########################################################################
##########################################################################
#########################################################################
# This snippet can be adapted to get rid of unwanted (NaN)values if required.
'''
for i in range(0, npoints_theta):

    if CoBeam_test[i] > 0:
        CoBeam_test_dB[i]=10*np.log10(CoBeam_test[i]) 
    else:
         CoBeam_test_dB[i] = -10000


'''



###############################

print("THIS IS THE END !!! (my friend)")    
 