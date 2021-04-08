# HGP2 - A (hopefully) fast, simple implementation of Head-Gordon, Pople's [refs]
# ERI recurrance relations.
#
# We're going to break this into several steps:
#
# 1. Primitive shell generation [ab,cd]
# 
# Given a shell-level description of a basis set ([ss,ss], [sp,ss], [dd,df], 
# including l [ll,ss]) generate all ERIs for the primitive basis functions in that shell.
# 
# 2. Contraction to (ab,cd)
# 
# 3. Integral Array Generation
# 
# Populate an integral record with the relevant terms for use in 
# an electronic structure theory code.
# 
# Gill's work on PRISM suggests [refs] being more flexible about when the basis function
# contraction is performed can reduce the operations count, but this will be simple and 
# likely fast.
# 
# There are many other recurrence relations to consider (MD [ref], LRL [ref], 
# OS [ref], Rys [refs]) but this should be a template for those others.
# 


# 0. Warm up
#    To see how well this works in practice, it might be useful to generate a few
#    simple cases and run them through steps 1-3 to see what's wrong with the plan.
#    Therefore, here are a few simple warm up exercises:
#
#    A. (00,00) generation aka (ss,ss)
#    B. (11,11) generation aka (pp,pp)
#    C. Integral array generation for h2o/sto-3g, which should be do-able with the above.

# 1. Primitive shell generation [ab,cd]
#
#    A. [0]m generation
#    B. [a+b,c+d] generation (VRR)
#    C. [ab,cd] generation (HRR)
#

# 2. Contraction to (ab,cd)
# 
# 3. Integral Array Generation