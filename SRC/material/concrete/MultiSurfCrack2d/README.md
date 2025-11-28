# multi-surf-crack2D
Plasticity implementation of interface material model

## Brief Description
This material model, based in plasticity theory, has been developed to describe the cyclic response of cracks in reinforced concrete. 

Notable features of the model include:
 - multiple-yield-surface formulation
 - crack-width-based yield surfaces
 - configuration dependent behavior, i.e., distinction between aggregate interlock and "free slip"
 - valid for general cracks, including flexural cracks, shear cracks, mixed-mode cracks, and crushing
 
The material model is intended for use with a zero-length node-to-node interface element, and has been developed within the OpenSees framework

## Instructions for Compiling within OpenSees
We plan to submit our model to the official OpenSees repository, but until that time, the code can be compiled into local versions of OpenSees. The C++ source files need to be included into the OpenSees source code. Detailed instructions can be found at:

[text](https://github.com/RESSLab-Team/OpenSees-Instructions)
