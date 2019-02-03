# VBParafac2
Implementation of the probabilistic Parafac2 model described in *INSERT PAPER TITLE HERE*.

# Usage
*DESCRIPTION OF USING THE CODE HERE*


# Dependencies

## OVPCA
from the Variational Bayes section here:
http://staff.utia.cz/smidl/index.php?module=pagemaster&PAGE_user_op=view_page&PAGE_id=5&MMN_position=4:4


## MATLAB File Exchange
- logdet
- multinv
- mtimesx
- parafac2
- addaxis6
- linspecer
- tight\_subplot

# Setup

## mtimesx
should be compiled with the mex function.

try `mex -v -largeArrayDims mtimesx.c -lmwblas`

if that doesn't work then choose the supported version of gcc (if the default isn't supported a warning will printed) and set the `blas_lib` variable to `PATH_TO_MATLAB/bin/glnxa64/libmwblas.so`

then run

`mex('GCC="/usr/bin/gcc-4.9"','-DEFINEUNIX','mtimesx.c',blas_lib)`
