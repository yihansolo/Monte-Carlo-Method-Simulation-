#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 13:20:53 2018

@author: yihangao
"""

#Import modules
import random;
import numpy;
import math;

#Global Variables
E_Crit = 45; #Cirital energy in MeV
T_P = 0.9525; #Thickness of the plate in cm 
R_Det = 7.62; #Radius of detector in cm
X_0 = 8.9; #Radiation length in cm


#Function cacluates length particle will travel given an energy in MeV
def StoppingLen(E_Val):
    L_Stop = X_0*math.log(1 + E_Val/E_Crit);
    return L_Stop; #Stoping length in cm

#Number of sparks
def N_Spark(L,Theta,Z_0):
    Num_Spark = 1 + math.floor(L*math.cos(Theta) - Z_0/(T_P));
    return Num_Spark;
