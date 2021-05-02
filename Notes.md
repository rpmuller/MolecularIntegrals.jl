# Notes on changes to MolecularIntegrals.jl
## 2021-05-02
Attempted to pass in global arrays explicitly to local scope, i.e.:
```
function vrr(amax,cmax, aexpn,bexpn,cexpn,dexpn, A,B,C,D
    # Pass in global arrays explicitly to local scope
    ,shift_index=shift_index,shift_direction=shift_direction
    )
```
Didn't have any impact on timings (global ethane/6-31G/cc-pvdz=1.11/16.1, local 1.12,16.6)


