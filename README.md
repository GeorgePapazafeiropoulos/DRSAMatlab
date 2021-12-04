# Dynamic-Response-Spectrum-Analysis
Dynamic response spectrum analysis for estimation of structural response to dynamic loads (earthquakes, vibrations, shocks, etc.)


Matlab code for the application of the dynamic response spectrum analysis (DRSA) of structures is presented. This procedure utilizes the following steps:

    Calculation of the stiffness matrix of the structure for the free degrees of freedom
    Calculation of the (lumped) mass matrix of the structure
    Definition of the degrees of freedom which convey the dynamic loading in the structure
    Definition of the response spectra that will be used for the calculation of the structural response according to the DRSA. Five types of spectra can be used: pseudoacceleration, pseudovelocity, displacement, velocity and acceleration response spectra. These spectra can be either design spectra according to the various norms and guidelines, or actual earthquake spectra.
    Execution of the dynamic response spectrum analysis using the function DRSA.m.
    Calculation of various structure-dependent quantities such as base shear and base overturning moment in shear buildings.
    Application of a combination rule for the estimation of the maximum structural response (displacement, velocity, acceleration, interstorey drift, base shear, base moment, etc.). Three combination rules are available in the present code: absolute sum (implemented by the function ABSSUM), square root of the sum of squares (implemented by the function SRSS) and complete quadratic combination (implemented by the function CQC). 

The function DRSA.m utilizes the following steps:

    Calculate the eigenfrequencies and eigenmodes of the structure
    Calculate the generalized masses
    Calculate the participation coefficients for each eigenmode
    Calculate the response spectrum values for each eigenperiod of the structure by performing linear interpolation based on the response spectrum values that are input to the function
    Return the contribution of each eigenmode to the total response.

The present code is accompanied by 5 examples in which its application is presented. These examples are taken from various standard textbooks or other material. The results of the examples are verified by the results of the application of the present code.
Additional plots are added to the examples to graphically illustrate to the readers the whole dynamic response spectrum analysis procedure. 
The author is open to any suggestions or recommendations that the users may have.
