# TODO and notes

## Start a docs folder
- [X] Move general notes from HGP into docs
- Learn standard way of doing Julia docs
- Links to other integral projects, pyquante versions, and
    existing julia qchem packages.
- Expand docstrings to make them more useful.
- Put badge-link to docs in README.md file
- Add ref to Rys work
- Add doi numbers to refs

## Write hrr_dict and hrr_array versions
- [X] These are the old version of hrr1 and hrr5, respectively
    - [X] Will have to recover the old version of hrr1 from before 4/22.
- Should be able to call *either* vrr_array or vrr_wide
    - E.g. :use_widearray symbol/keyword option
- Comment code to indicate that we're supporting multiple
    interfaces for convenience, and that the hrr_array is
    faster, but the hrr_dict can only return the requisite
    terms.
- Prune out unnecessary results from hrr_dict (when != ashell, != bshell, etc.)
- [X] Rewrite timing routines to call new names
- Test/time/commit/push

## Write fuzz-like test functions in runtests.jl
- vrr_test(ashell,cshell) code:
    - Test all elements of vrr_array and vrr_widearray against ERI.jl code
    - Random A,C, aex, cex? Or maybe fallback if A,C,aex,cex not specified
    - Can use this as an excuse to write functions to generate random PGBFs
        in some distance/exp/powers
- hrr_test(ashell,bshell,cshell,dshell)
    - Test all elements of hrr_array and hrr_dict against ERI.jl code
    - Random A,B,C,D??? Or maybe fallback if A,B,C,D not specified
- Add comments to ERI code to the effect that it is slow code that is very
    likely to be correct that is kept around for reference and testing.
- Merge Low level OneInts tests with OneInts
- ERI tests should be against standard (large R) or pyquante values.
- vrr/hrr tests should be against ERI versions

## Move contracted routines to CHGP.jl
- Test chrr
- Write out a standard integral record using this code
- Write integral records to .jld files using JLD2 HDF5 format


## Release version 0.1.0
- Register julia package wherever I'm supposed to do this
- Post to discord


## Get SP shells working
- May have to redo indexing of shell_indices to make this work
- SP shells could be one of the things that is easier to do in the 
    generated code (below)

## Write HGPgenerate.jl script
- Julia code that generates more julia code
- Will generate a HGPgen.jl file that contains generated hrr_abcd and vrr_ab code
- Fall back to hrr_array and vrr_array versions
- Which directory do support scripts like this live in? tools?
- Write test routines to test **all** of HGPgen.jl.
- Expand timing routines to include HGPgen routines.
- Are there efficiencies I can use to compile HGPgen.jl code directly
    once I generate it?


## Include basis sets in g94 format
- Don't need to include every one in BSE
- Check against files in Data.jl
- Can replace SP with S,P if SP shells aren't working yet.
- Is this in a /data directory? How do I support it?


## Bump version to 0.2




# Someday changes
- Make bf I,J,K and amplitude() calls arrays [I,J,K]
- Decide on a better m2ao: are the arguments arrays or tuples. Are they both?
- Overhaul m2ao (support tuples & arrays as keys??)

## Derivatives??
- Look at what derivative code julia routines support.

## Reach out to julia package maintainers to use MolecularIntegrals?
- GAMESS team, psi4 team, etc.

## Implement other molecule methods:
- nocc, nclosed, nopen, nup, ndown, stoich, center!



# Completed tasks kept for records

## Comment HGPold routines
- [X] Indicate that they're kept for reference, but no longer
    included in codebase or tested
- [X] Remove hrr2/vrr2 from timing.jl

## Write vrr_array and vrr_widearray versions
- [X] These are vrr5 and vrr1, respectively
- [X] Comment the code to indicate that we're supporting multiple
      interfaces for convenience, and that that the speed is 
      roughly equivalent
- [X] Rewrite timing/testing routines to call new names
- [X] Test/time/commit/push