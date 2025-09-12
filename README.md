# Cross-bridge model-based quantification of muscle metabolite alterations leading to fatigue during all-out knee extension exercise 



## Abstract

Intense physical exercise is associated with high energy demands and muscle metabolite changes that affect force generation, leading to muscle fatigue. Although these changes are well characterized in humans, their contribution to muscle fatigue is not clearly understood. Furthermore, we lack methodologies for a systems-level exploration of these changes that occur during intense exercise to understand the mechanisms behind muscle fatigue development. Therefore, in this study, we developed an updated human skeletal muscle model and included new proton-binding mechanisms and accounted for the kinetics of cross-bridge cycling and its associated metabolite-mediated inhibition as well as key metabolic processes associated with ATP production. We contextualized and parameterized the updated model using human muscle activation data collected during an all-out knee extension exercise as input to simulate the alterations in force generation and the temporal evolution of muscle metabolite levels (Pi, H+, and PCr) as outputs. Our model predicted that the experimentally observed progressive decline in muscle activation during exercise had a minimal effect on force generation and that the accumulation of Pi and H+ inhibited force generation both individually and synergistically. Furthermore, we found that force generation was more sensitive to H+ than Pi during an all-out knee extension exercise, with Pi impacting force generation by increasing actin-myosin detachment and H+ by preventing the formation of a strongly bound cross-bridge state. Our computational analysis was able to identify the mechanisms behind muscle force generation and quantified the potential roles of muscle activation, Pi, and H+ changes in fatigue development during an all-out knee extension exercise. 

# Contents

This repository contains the data, scripts, and models used in this study. They are organized into two folders 1) raw\_data and 2) codes.

* Contents of folder /raw\_data are as follows:

  * Power, EMG, Phosphocreatine concentration, Pi concentration , ADP concentration and pH dataset used for parameterization are in the excel files power\_for\_fitting\_DPF\_2.xlsx, Emg\_for\_fitting\_DPF.xlsx, Pcr\_for\_fitting\_DPF.xlsx, pi\_for\_fitting\_DPF.xlsx, ADP\_for\_fitting\_DPF.xlsx, and pH\_for\_fitting\_DPF.xlsx, respectively.
  * Resting state concentrations and the sarcomere shortening velocity data used for parameterization are provided in the excel files Initial\_state.xlsx and dsdt\_for\_fitting\_DPF\_2.xlsx, respectively.
  * The data set used for validation is provided in the folder val\_dataset

* Contents of folder /codes are as follows:

  * This folder contains the MATLAB codes (MATLAB/R2022a) that represent the musculoskeletal model and the scripts used to simulate the model results shown in this manuscript
  * Model\_XB\_human\_QC.m encodes the 4-state crossbridge model and the kinetic model discussed in the manuscript and is used to simulate the force and metabolite dynamics
  * params.xlsx in folder /params contains the model parameters estimated using our parameterization routine and subsequently used to simulate the model
  * figure\_4\_subplots.m simulates the model in Model\_XB\_human\_QC.m  to generate the subplots of Figure 4
  * figure\_5\_subplots.m simulates the model in Model\_XB\_human\_QC.m  to generate the subplots of Figure 5
  * figure\_6\_a.mlx, figure\_6\_b.mlx, figure\_6\_c.mlx and figure\_6\_d.mlx simulates the model in Model\_XB\_human\_QC.m to generate the subplot A, B, C and D in Figure 6
  * figure\_8\_b.mlx, figure\_8\_c.mlx, figure\_8\_d.mlx and figure\_8\_e.mlx simulates the model in Model\_XB\_human\_QC\_SI.m to generate the subplot B, C, D and E in Figure 8. Model\_XB\_human\_QC\_SI.m is the model that incorporates the alternate proton inhibition hypothesis depicted in Figure 8A.
  * figure\_s2\_sublots.m simulates the model in Model\_XB\_human\_QC\_SI.m  to generate the subplots of Figure S2.

# Set up

* The codes were written and tested in MATLAB/R2022a and they can be directly executed from MATLAB command prompt or from the editor

Example: To run the code figure\_4\_subplots.m, directly type 'figure\_4\_subplots' in the MATLAB command prompt or open the editor and run it by pressing F5 or clicking the 'Run' icon in the toolbar of MATLAB Editor.





