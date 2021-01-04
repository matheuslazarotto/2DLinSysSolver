#!/bin/bash

######################################
#       2D Linear System Solver      #
###################################### 
# This code runs and plots the numeri-
# cal solution for a given 2D linear 
# dynamical system.
#
# Instructions:
# => Set matrix coeffientes in varia-
# bles [A11, A12, A21, A22].
# 
# => Set phase-space interval window 
# to show (initial conditions are ta-
# ken inside this window) in variables 
# [x_min, x_max, y_min, y_max]
# 
# => Set number of initial conditions 
# along each axis [Nx, Ny]
# 
# => In case any change to 'main.cpp' 
# code is made, the flag [compile] must 
# be selected as 'true'.
#
# Implemented by: M. Lazarotto (30/10/2020)

######################################
#              Settings              #
######################################
# Matrix [A] coefficients
A11="1.0"                            # 
A12="1.0"                            # A = (A11  A12)
A21="-1.0"                           #     (A21  A22)
A22="-1.0"                           # 
# Phase-space (x,y) limits
x_min="-1.5"                         # Min x val
x_max="1.5"                          # Max x val
y_min="-1.5"                         # Min y val
y_max="1.5"                          # Max y val
# Initial conditions
Nx=50                                # Points along x
Ny=1                                 # Points along y
# Integration time 
t_run="50.0"                         # Total time to evolve
# Compilation flag
compile=true                         # Select to compile
######################################
#                                    #
######################################

if [ "$compile" = true ] || [ ! -e "./2dlinsyssolver" ]; then
    if [ -e "./2dlinsyssolver" ]; then
        make clean
    fi
    make
fi

if [ -e "./2dlinsyssolver" ]; then
    clear
    
    # Run dynamics
    ./2dlinsyssolver $A11 $A12 $A21 $A22\
                     $x_min $x_max $y_min $y_max\
                     $t_run $Nx $Ny

    # Plot result
    python3.6 python/plot.py -x_min "$x_min" -x_max "$x_max"\
                             -y_min "$y_min" -y_max "$y_max"

    date
fi