#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 13:20:53 2018

@author: yihangao
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 14:44:24 2018

@author: jerem
"""


#Import modules
#Import all needed modules
import random;
import numpy;
import math;
import array;
import matplotlib.pyplot as plt

#Global Variables
Pi = math.pi;
E_Crit = 45; #Cirital energy in MeV
T_P = 0.9525; #Thickness of the plate in cm 
T_G = 0.635; #Thickness of gaps in cm
R_Det = 7.62; #Radius of detector in cm
R_Fid = R_Det/2; #Maxium radius of electrons in cm
X_0 = 8.9; #Radiation length in cm
Theta_Max = Pi/6; #Maxium Polar angle
N_Exp = 43; #Number of muon decays observed experimentally
N_Run = 10000; #Total Number of runs per simulation
N_exp_array=[4,7,13,16,4] 
S=[]
 
 #Function cacluates length particle will travel given an energy in MeV
def StoppingLen(E_Val):
     L_Stop = X_0*math.log(1 + E_Val/E_Crit);
     return L_Stop; #Stoping length in cm
 
 #Number of sparks
def N_Spark(L,Theta,Z_0):
     Num_Spark = 1 + math.floor((L*math.cos(Theta) - Z_0)/(T_P));
     return Num_Spark;
 
 #Escape length function (calculates escape length for given theta, r and phi)
def Escape_L(r,phi,theta):
    L_escp = (-r*math.cos(phi) + math.sqrt(r**2 * (math.cos(phi))**2 + (R_Det**2 - r**2)))/(math.sin(theta))*(T_P/(T_P + T_G));
    return L_escp;
 #Define N_Sim function
def NSim(rad,E_e,m_mu):
    N_Val = rad*(m_mu*E_e)**2 *(3 - (4*E_e)/(m_mu))
    return N_Val #Returns un normalized N_Sim Value
def res_sum(AssociatedDecay):
    sum=0 
    for k1 in range(2,len(AssociatedDecay)-1):
        res_sq=(N_exp_array[k1-2]-AssociatedDecay[k1])**2
        sum=sum+res_sq
    return res_sq
    
    
 #Generate array of muon masses in MeV
LowMuonMass = 80; #Lowest Muon mass
HighMuonMass = 110; #Highest muon 
MassRuns = (HighMuonMass - LowMuonMass) + 1; #Number of dimassfferent masses to be tested  
MuonMass_Array = numpy.linspace(LowMuonMass,HighMuonMass,MassRuns); #Muon mass array increments by 1 MeV (All masses are in MeV)
  #These matricies will keep track of the data from the simulations, each column corresponds to a different muon mass and each row corresponds to 
 #run number
Spark_Matrix = numpy.zeros((N_Run,MassRuns)); #Initialize an empty matrix corresponding to number of sparks
Energy_Matrix = numpy.zeros((N_Run,MassRuns)); #Initialize an empty matrix corresponding to electron energies
Radius_Matrix = numpy.zeros((N_Run,MassRuns)); #Initialize an empty matrix corresponding to electron radius
MuonDec_Matrix = numpy.zeros((N_Run,MassRuns)); #Initialize an empty matrix corresponding to unnormalized N_Sim
 #Debug_Matrix = numpy.zeros((N_Run,5)); #For debugging code
  #Simulation loop
for i1 in range(0,MassRuns): #Loop over every muon mass
    for i2 in range(0,N_Run): #Run the simulation with a given muon of mass N times
         #Generate Electron Energy
         Electron_En = random.uniform(0,0.5*MuonMass_Array[i1]); #Energy is in MeV
         Energy_Matrix[i2,i1] = Electron_En; #Add electron energy to matrix
         
         #Generate electron radius, polar angle, azimuthal angle and thickness
         r_Val = random.uniform(0,R_Fid);
         Radius_Matrix[i2,i1] = r_Val; #Add radius to matrix
         Theta_Val = random.uniform(0,Theta_Max); #Generate polar angle
         Phi_Val = random.uniform(0,2*Pi); #Generate azimuthal angle
         Z_Val = random.uniform(0,T_P); #Generate a thickness in cm
         
         #Calculate length and escape length
         Len_Val = StoppingLen(Electron_En);
         Escp_Len = Escape_L(r_Val,Phi_Val,Theta_Val);
        
        #Obtain unormalized N_Sim Value
         MuonNum = NSim(r_Val,Electron_En,MuonMass_Array[i1]);
         MuonDec_Matrix[i2,i1] = MuonNum;
        
        #Compare escape length and Length
         if (Len_Val > Escp_Len):
            Len_Val = Escp_Len; #Set length to escape length
         
         
        #Calculate the number of sparks produced and add to data matrix
         Spark_Num = N_Spark(Len_Val,Theta_Val,Z_Val);
         Spark_Matrix[i2,i1] = Spark_Num;
        
        
    
    #Normalize NSim Matrix row before incrementing the mass
    EventSum = 0; #Sum of unormalized muon detection events
    for i3 in range(0,N_Run):
        EventSum = EventSum + MuonDec_Matrix[i3,i1];
    #Get the normalization constant
    D = N_Exp/(EventSum);
    #Scale the Muon Decay column by the calculated normalization factor
    MuonDec_Matrix[:,i1] = D*MuonDec_Matrix[:,i1];    
           
        
##Data Analysis Section##
#For each number of sparks determine how many muon decays occured
#for j0 in range(0,MassRuns)
 #Which mass run is being considered
for j0 in range(0,MassRuns):
    SparkArray = [1,2,3,4,5,6,7,8];
    AssociatedDecay = array.array('d');
    for j1 in range(0,8): #Loop over all possible numbers of sparks
        IndxArray = array.array('d'); #Location of all j1 sparks in the spark matrix
        #Find the indicies of all n_Spark = j1
        for j2 in range(0,N_Run):
            if (j1 == Spark_Matrix[j2,j0]):
                IndxArray.append(j2);
        
        #Sum over all the indices in the muon decay array to determine the number of decays with n_spark
        Muon_Sum = 0;
        for j3 in range(0,len(IndxArray)):
            IndVal = int(IndxArray[j3])
            Muon_Sum = Muon_Sum + MuonDec_Matrix[IndVal,j0];
        AssociatedDecay.append(Muon_Sum); #Append the number of muon decays
        print(Muon_Sum)
    #Calculate squared residuals 
    sum=res_sum(AssociatedDecay)
    S.append(sum)


#Experimental Data 
ExpSpark = [3, 4, 5, 6, 7]                  # number of sparks
ExpDecayRate = [4, 7, 13, 16, 4]            # number of decays
ExpErrorBars = [2, 2.65, 3.61, 4, 2]        # length of error bars

#Plot the number of muon decays vs number of sparks
plt.figure(1)
plt.bar(SparkArray,AssociatedDecay, color = (0.8,0.2,0.2), label='Simulation', zorder=1)
plt.title("")
plt.xlabel("Number of Sparks");
plt.ylabel("Number of Muon Decays");
plt.scatter(ExpSpark,ExpDecayRate, color = (0,0,0), label='Experimental Data', zorder=3)
plt.errorbar(ExpSpark,ExpDecayRate, yerr=ExpErrorBars, color = (0,0,0), linestyle='none', zorder=2)
plt.legend(loc='upper left')

plt.figure(2)
x=[]
for k2 in range(80,111):
    x.append(k2)
plt.scatter(x,S)

plt.show()
