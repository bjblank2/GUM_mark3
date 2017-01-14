# GUM_mark3

This python-based platform predicts the magnetocaloric performance of 
magnetic Heusler shape memory alloys.  

Approach: Shape memory alloy properties are described by an 
integrated Blume-Emery-Griffiths description of the martensitic phase 
transition and a cluster-expansion.  The former captures the 
entropy-driven martensitic phase transformation, while the latter 
captures chemical interaction and magnetic exchange coefficients. 

There are two parts to the platform:

(1) Model fitting: Fit the effective model parameters to the results of 
ground state DFT calculations, using cross-validation to optimize the set 
of clusters used in the expansion.

(2) Monte Carlo: Simulation tools to predict temperature dependent properties 
using the fitted model parameters. 

Main Developer: Brian J. Blankenau
Contributors: Chendi Lin, Elif Ertekin.
