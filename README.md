# Tracing information propagation in a fish school with a conditioned escape response

This repository contains the code to reproduce the results from the manuscript: Jadhav, V., 
Aimon, C., Lamshana, F., Trendafilov, D., Escobedo, R., Sire, C., Theraulaz, G., \& 
Guttal, V. (2026). Tracing information propagation in fish schools with a conditioned 
escape response. 

The codes are tested to run on MATLAB Version: 25.2.0.3177638 (R2025b). Ensure that the `Signal processing` and `Econometrics` toolboxes are installed. 

## Data

Tracked fish trajectory data used in Figure 1e are available in the `/main_text/` folder. Tracks 
of each fish are labelled as `Clip0036_1C_4N_E4_4_bee1.csv`, `Clip0036_1C_4N_E4_4_bee2.csv`, $\dots$, 
`Clip0036_1C_4N_E4_4_bee5.csv`. The first, second, and third columns in the `.csv` files 
represent the frame number (at 100 frames per second), the position along the x-axis (cm), 
and the position along the y-axis (cm), respectively.

Tracked fish trajectories from all 53 trials are in `/main_text/all_bd_data.mat`. 
`all_bd_data.mat` contains the position (`pos_t_ex_x`), speed (`speed_t_ex_x`), 
velocity (`vel_t_ex_x`), and the frame at which the green light was turned on (`frame_light_on_ex_x`) 
for all experiments, where `x` is the trial number from 1 to 53.

The data file `all_bd_data.mat` must be in the `/main_text/` folder to reproduce all figures 
in the main text and electronic supplementary material text. This file is already included in 
the folder.

## Code for figures 1e and 2a, b, c

**Trajectories of individual fish, speed, polarisation, and dispersion from a representative 
trial of the collective escape response:** run `/main_text/figure_1e_2ace.m` to generate figure 1e 
and figure 2a–c.

## Code for figure 2d–f, figure 3, and figure 4a–b

Run `main_text/col_beh_ana.m`: This code calculates the average polarisation, dispersion, and 
speed as a function of normalised time *t*. It also calculates the probability density functions 
of polarisation, dispersion, and speed across the initial, escape, and relaxation phases. Finally, 
it calculates the normalised crossing time since the green light turned on and the 
normalised time since the previous fish crossed the hurdle, both as a function of crossing rank.

_Note:_ Ensure `getCor_scalar.m` is in the same folder. 

## Code for figure 4c

**Leader-follower interaction network during collective escape:** Run `main_text/fig_4c.m` to 
generate the leader-follower interaction network for the experimental data.

## Code for results from the model of collective escape dynamics

We provide the required `.csv` files containing the model results in the `/main_text` 
directory to reproduce all the figures in the main text. However, all the data can be 
generated using the code provided in the `/model` folder.

1. **Simulate the model:** Run `/model/sim_data.m` to simulate the model. 
2. **Calculate group properties:** Run `model/grp_properties.m` to generate figure 2g–h, 
figure 3, figure 4a–b, and electronic supplementary material figures S9 and figure S10. 
3. **Leader-follower interaction network:** Run `/model/pert_test_exp.m` to construct the 
leader-follower network described in Section 3b of the main text. This script 
generates a `.csv` file containing the leader-follower data. Run `/model/interaction_net_ana.m` 
to generate figure 4c. 
4. **Simulation video:** Run `/model/simulate.m` to visualise the collective escape dynamics 
obtained from the model.


## Code for results from the null model of collective escape dynamics

1. **Simulate the model:** Run `/null_model/sim_data_mp.m` to simulate the null model described 
in the electronic supplementary material section S1.
2. **Calculate group properties:** Run `/null_model/grp_properties.m` to generate electronic 
supplementary material figure S3. This code also generates the mean group speed for the null model (figure S6), 
which can be compared with figure 2g of the collective escape model. 
3. **Leader-follower interaction network:** Run `/null_model/pert_test_exp.m` to construct the 
leader-follower network described in Section 3b of the main text. This script 
generates a `.csv` file containing the leader-follower data. Run `/model/interaction_net_ana.m` 
to generate the network. In the null model, no leader-follower pairs are found because there 
are no interactions between agents.  
4. **Simulation video:** Run `/null_model/simulate.m` to visualise the group dynamics obtained 
from the null model.

## Supplementary figures

1. `/main_text/col_beh_ana.m` generates electronic supplementary material figure S1, S5
2. `/sm/fig_s4a.m` generates electronic supplementary material figure S4a
3. `/sm/fig_s4bc.m` generates electronic supplementary material figure S4b-c. 



