# FOXcombustionModels
This is a library containing alternative combustion models for OpenFOAM. These models are primarily intended for LES, but may work for RANS as well.
The library was made for OpenFOAM v2112, but you should be able to adapt it to any other version without too much effort.

Included models:


**foxPaSR**

This model attempts to scale reactions separately according to the PaSR approach, with individual time scales for all reactions. The turbulent mixing time scale is a harmonic average between the Kolmogorov time scale and the time scale of sub-grid velocity stretch. The model is described in detail in Åkerblom 2023 (10.17196/OS_CFD#YEAR_2022) NOTE: CURRENTLY ONLY WORKS FOR SINGLE-STEP MECHANISMS. I hope to remedy this in the future.

**renPaSR**

This model is similar to the default PaSR model, but computes the time scales differently. The chemical time scale is based on the shortest species residence time, following Ren & Golding 2011 (10.1016/j.combustflame.2011.02.018). The turbulent mixing time scale is a harmonic average between the Kolmogorov time scale and the time scale of sub-grid velocity stretch.
