# TODO and notes

## Timing and benchmarks
- [ ] Compare to JERI.jl benchmarks
- [ ] Would be nice to be within a factor of 10 by release 0.1.0

## Contracted routines
- [X] Determine whether the normalization constants are the same for all m-values corresponding to an L-value. They're NOT.
- [X] Write out a standard integral record using this code
- [X] Time for ethane
- [ ] Figure out right place to add normalization consts
- [ ] Working test chrr

## More documentation
- Host a webpage for MolecularIntegrals.jl
- Put badge-link to docs in README.md file

## Release version 0.1.0
- Register julia package wherever I'm supposed to do this
- Post to discourse server

## Write HGPgenerate.jl script
- Julia code that generates more julia code
- Will generate a HGPgen.jl file that contains generated hrr_abcd and vrr_ab code
- Fall back to hrr and vrr versions
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

## Write fuzz-like test functions in runtests.jl
- vrr_test(ashell,cshell) code:
    - Test all elements of vrr against ERI.jl code
    - Random A,C, aex, cex? Or maybe fallback if A,C,aex,cex not specified
    - Can use this as an excuse to write functions to generate random PGBFs
        in some distance/exp/powers
- hrr_test(ashell,bshell,cshell,dshell)
    - Test all elements of hrr and hrr_dict against ERI.jl code
    - Random A,B,C,D??? Or maybe fallback if A,B,C,D not specified
- Add comments to ERI code to the effect that it is slow code that is very
    likely to be correct that is kept around for reference and testing.
- Merge Low level OneInts tests with OneInts
- ERI tests should be against standard (large R) or pyquante values.
- vrr/hrr tests should be against ERI versions

## Bump version to 0.2

## Get SP shells working
- May have to redo indexing of shell_indices to make this work
- SP shells could be one of the things that is easier to do in the 
    generated code (below)

# Someday changes

## Derivatives??
- Look at what derivative code julia routines support.

## Reach out to julia package maintainers to use MolecularIntegrals?
- GAMESS team, psi4 team, etc.

## Implement other molecule methods:
- nocc, nclosed, nopen, nup, ndown, stoich, center!
