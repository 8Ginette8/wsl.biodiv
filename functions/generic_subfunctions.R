# ###########################################################################
# Generic second order functions
#
# $Date: 2018-05-08 V2

# Author: Philipp Brun, philipp.brun@wsl.ch
# Dynamic Macroecology Group
# Swiss Federal Research Institute WSL
# 

# ###########################################################################

### =========================================================================
### define hde function
### =========================================================================
# hide annoying prints

hde=function(x){
  invisible(capture.output(x))
}