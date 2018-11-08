#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 13:20:53 2018

@author: yihangao
"""

#Import modules
from  __future__ import division
import random;
import numpy;
import math;

#Global Variables
Pi = math.pi;
E_Crit = 45; #Cirital energy in MeV
T_P = 0.9525; #Thickness of the plate in cm 
T_G = 0.635; #Thickness of gaps in cm
R_Det = 7.62; #Radius of detector in cm
R_Fid = R_Det/2; #Maxium radius of electrons in cm
X_0 = 8.9; #Radiation length in cm
Theta_Max = Pi/6; #Maxium Polar angle


N_Run = 100; #Total Number of runs persimulation


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
    L_escp = -r*math.cos(phi) + math.sqrt(R_Det**2 - r**2)/(math.sin(theta)*(T_P/(T_P + T_G)));
    return L_escp;
    

#Generate array of muon masses in MeV
LowMuonMass = 90; #Lowest Muon mass
HighMuonMass = 105; #Highest muon mass
MassRuns = (HighMuonMass - LowMuonMass) + 1; #Number of different masses to be tested  
MuonMass_Array = numpy.linspace(LowMuonMass,HighMuonMass,MassRuns); #Muon mass array increments by 1 MeV (All masses are in MeV)

#These matricies will keep track of the data from the simulations, each column corresponds to a different muon mass and each row corresponds to 
#run number
Spark_Matrix = numpy.zeros((N_Run,MassRuns)); #Initialize an empty matrix corresponding to number of sparks
Energy_Matrix = numpy.zeros((N_Run,MassRuns)); #Initialize an empty matrix corresponding to electron energies
Radius_Matrix = numpy.zeros((N_Run,MassRuns)); #Initialize an empty matrix corresponding to electron radius

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
        Escp_Len = Escape_L(r_Val,Theta_Val,Phi_Val);
        
        #Compare escape length and Length
        if (Len_Val > Escp_Len):
            Len_Val = Escp_Len; #Set length to escape length
        
        #Calculate the number of sparks produced and add to data matrix
        Spark_Num = N_Spark(Len_Val,Theta_Val,Z_Val);
        Spark_Matrix[i2,i1] = Spark_Num;
