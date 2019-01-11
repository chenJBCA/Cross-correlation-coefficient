#Instructions:

#Purpose: To compute the cross-spectrum of two maps and return their cross-correlation coefficient

#Defination: Cross-correlation coefficient = cross-spectrum/sqrt(auto-spectrum1 X auto-spectrum2)

#Requirement: IDL; Polspice 3.03; Healpix 3.11;

#Input files: 2 maps to be cross-correlated; 1 weight map to be applied to the cross-correlation.

#Power spectrum estimator: Primarily with polspice; Optionally using ianafast.

#Output polspice fits files: map1 auto-spectrum; map2 auto-spectrum; map1Xmap2 cross-spectrum; map1Xmap2 cross-correlation coefficient.

#Optional anafast output: map1 auto-spectrum; map2 auto-spectrum; map1Xmap2 cross-spectrum; map1Xmap2 cross-correlation coefficient.

#To run: 1.Type 'idl' in your terminal; 2.Type '.run szcor.pro'; 3. Type 'szcor'.
