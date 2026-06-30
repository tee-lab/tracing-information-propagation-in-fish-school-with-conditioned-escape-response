# Tracing information propagation in a fish school with a conditioned escape response

This repository contains the code for the manuscript:

The codes are tested to run on MATLAB Version: 25.2.0.3177638 (R2025b).

## Raw data

Tracked trajectory data used in figure 1e are available in the folder `/main_text/`. Tracks of 
each fish are labelled as `Clip0036_1C_4N_E4_4_bee1.csv`, `Clip0036_1C_4N_E4_4_bee2.csv`, ..., 
`Clip0036_1C_4N_E4_4_bee5.csv`. The first, second, and third columns in the `.csv` files are 
Frame number, position along the x (cm) and y(cm) axis, respectively. 

Tracked data from all the 53 trials are in `/main_text/all_bd_data.mat`. `all_bd_data.mat` contains
position (`pos_t_ex_x`), speed (`speed_t_ex_x`), velocity (`vel_t_ex_x`), and
frame at which green light was turned on (`frame_light_on_ex_x`) for all the experiments, where 
`x` is the trial number between 1 to 53. 

Data file `all_bd_data.mat` should be in `\main_text` folder to reproduce all the figures 
in the main text and supplementary text. We have already added the file in the folder.

## Code for figure 1e and 2a,b,c

**Trajectories of individual fish, speed, polarisation and dispersion from the representative 
trial of collective escape response:** Run `\main_text\figure_1e_2ace.m` to generate figure 1e and
figure 2a-c.

## Code for figure 2 d-f, figure 3, and figure 4 a-b

Run `main_text\col_beh_ana.m`: this code calculates average polarisation, dispersion and speed
as a function of normalised time *t*. It further calculates probability density functions of 
polarisation, dispersion and speed across initial, escape and relax phases. It finally calculates 
normalised crossing time since green light is turned and normalised time since previous fish crossed
the hurdle as a function of crossing rank. 

## Code for figure 4c

**Leader-follower interaction network during collective collective escape:** 
run `main_text\fig_4c.m` to generate the leader-follower interaction for experimental data.

## Codes for results from model of collective escape dynamics



