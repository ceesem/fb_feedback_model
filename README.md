Code for the rate models used in the mushroom body feedback paper.
The code should work with any remotely current version of MATLAB,
but requires export_fig and cbrewer, both of which can be found on
File Exchange. By Casey Schneider-Mizell (caseysm@gmail.com), 2018.

---
File description:

`all_feedback_models.m`

This set of scripts sets the neuronal, network connectivity, and stimulation
parameters for running the feedback model and plotting the results as seen in
the main text and supplemental figures of the paper.

`logistic_integration_general.m`

This file contains a generic rate model solver using ode45 and some associated
response shape functions.
