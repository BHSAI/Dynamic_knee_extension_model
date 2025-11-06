# Cross-bridge model-based quantification of muscle metabolite alterations leading to fatigue during all-out knee extension exercise



## Abstract

Intense physical exercise is associated with high energy demands and muscle metabolite changes that affect force generation, leading to muscle fatigue. Although these changes are well characterized in humans, their contribution to muscle fatigue is not clearly understood. Furthermore, we lack experimental methodologies for a systems-level exploration of these changes that occur during intense exercise to understand the mechanisms behind muscle fatigue development. In this study, we updated our previously developed human skeletal muscle model to include new proton-binding mechanisms and adapted it to study fatigue development during an intense all-out knee extension exercise. We contextualized and parameterized the updated model to simulate muscle force generation and muscle metabolite alterations, using muscle activation data obtained from human subjects performing an all-out knee extension exercise. Our model predictions showed that nullifying the observed decline in muscle activation during all-out exercise was not sufficient to stop fatigue development and suggested that other factors may play a role. We found that the accumulation of inorganic phosphate (Pi) and protons (H+), both individually and synergistically, were the main contributing factors at the cross-bridge level that inhibited force generation during all-out exercise. Furthermore, our model simulations showed that force generation was more sensitive to H+ than Pi during an all-out knee extension exercise, with elevated Pi levels promoting actin-myosin detachment and elevated H+ levels preventing the formation of strongly bound cross-bridge states. Our combined experimental/computational analysis also showed that both exercise routine and the type/group of muscles involved in muscle force generation determined the potential contributing factors responsible for fatigue development.



# Contents

This repository contains the data, scripts, and models used in this study. They are organized into two folders 1) raw\_data and 2) codes.

* Contents of folder /raw\_data are as follows:

  * Force, EMG, Phosphocreatine concentration, Pi concentration , ADP concentration and pH dataset used for parameterization are in the excel files force\_data\_broxterman.csv, iEMG\_data\_burnley.xlsx, PCr\_data\_broxterman.csv, Pi\_data\_broxterman.csv, ADP\_data\_broxterman.xlsx and proton\_data\_broxterman.xlsx respectively.
  * Resting state concentrations used for parameterization are provided in the excel files initial\_state.xlsx respectively.

* Contents of folder /codes are as follows:

  * This folder contains the MATLAB codes (MATLAB/R2022a) that represent the musculoskeletal model and the scripts used to simulate the model results shown in this manuscript
  * Model\_XB\_human\_QC.m encodes the 5-state crossbridge model and the kinetic model discussed in the manuscript and is used to simulate the force and metabolite dynamics
  * params.xlsx in folder /params contains the model parameters estimated using our parameterization routine and subsequently used to simulate the model
  * figure\_2\_subplots.m simulates the model in Model\_XB\_human\_QC.m  to generate the subplots of Figure 2
  * figure\_3\_subplots.m simulates the model in Model\_XB\_human\_QC.m  to generate the subplots of Figure 3
  * figure\_4\_a.m, figure\_4\_b.m, figure\_4\_c.m, figure\_4\_d.m, and figure\_4\_e.m simulates the model in Model\_XB\_human\_QC.m to generate the subplot a, b, c, d, and e of Figure 4
  * figure\_5\_a.m, figure\_5\_b.m, figure\_5\_c.m, figure\_5\_d.m, and figure\_5\_e.m simulates the model in Model\_XB\_human\_QC.m to generate the subplot a, b, c, d, and e of Figure 5

# Set up

* The codes were written and tested in MATLAB/R2022a and they can be directly executed from MATLAB command prompt or from the editor

Example: To run the code figure\_2\_subplots.m, directly type 'figure\_2\_subplots' in the MATLAB command prompt or open the editor and run it by pressing F5 or clicking the 'Run' icon in the toolbar of MATLAB Editor.

