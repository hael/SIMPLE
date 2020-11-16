#!/usr/bin/env julia
# we use strings for reading
ostr    = Array{String}(2)
ostr[1] = "   e1=0.00 e2=59.000 e3=6.0 x=5.5678 y=4.56677 state=1 movie=intg1.mrc "
ostr[2] = " e1=5.00 e2=178.000 e3=123.0 x=0.4578 y=2.56677 state=0 movie=intg2.mrc   "
# and array of dictonaries for storage/manipulation
os  = Array{typeof(Dict{String,Any}())}(2)
ori = Dict{String,Any}()
for i in 1:2
   ostr_split = split(strip(ostr[1]))
   for j in eachindex(ostr_split)
      keyval = split(ostr_split[j], "=")
      ori[keyval[1]] = keyval[2]
   end
   os[i] = ori
end
show(os)
