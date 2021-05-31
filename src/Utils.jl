using LinearAlgebra
using StaticArrays
export factorial2, dist2, ipairs, rpairs, spairs, iindex

@inline factorial2(n::Int64) = prod(n:-2:1) # double factorial !!
@inline dist2(dx,dy,dz) = dx*dx+dy*dy+dz*dz # TODO: use hypot()^2?
@inline dist2(dxyz) = dot(dxyz,dxyz)
@inline dist2(xyz1,xyz2) = dist2(xyz1-xyz2)

ipairs(n::Int64) = ((i, j) for i = 1:n for j = 1:i)
rpairs(n::Int64) = ((i,j) for i in 1:n for j in 1:n) # rectangular option to old pairs
spairs(n::Int64) = ((i, j) for i = 1:n for j = 1:(i-1)) # subdiagonal option to old pairs
 
triangle(i::Int64) = div(i*(i+1),2)
triangle(i::Int64,j::Int64) = i<j ? triangle(j-1)+i : triangle(i-1)+j
                        
@inline iiterator(n::Int64) = ((i,j,k,l) for (i,j) in ipairs(n) for (k,l) in ipairs(n) if triangle(i,j) <= triangle(k,l))

@inline iindex(i::Int64,j::Int64,k::Int64,l::Int64) = triangle(triangle(i,j),triangle(k,l))



