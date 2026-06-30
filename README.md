# Tracing information propagation in a fish school with a conditioned escape response

This repository contains the code for the manuscript:

The codes are tested to run on MATLAB Version: 25.2.0.3177638 (R2025b).

## Raw data

Tracked trajectory data used in figure 1e are available in the folder `/main_text/`. Tracks of 
each fish are labelled as `Clip0036_1C_4N_E4_4_bee1.csv`, `Clip0036_1C_4N_E4_4_bee2.csv`, ..., 
`Clip0036_1C_4N_E4_4_bee5.csv`. The first, second, and third columns in the `.csv` files are 
Frame number, position along the x (cm) and y(cm) axis, respectively. 

Tracked data from all the 53 trials are in `/main_text/all_bd_data.mat`. `all_bd_data.mat` contains
position, speed, velocity, and frame at which green light was turned on for all the experiments. 

## Code for figure 1e and 2a,b,c
