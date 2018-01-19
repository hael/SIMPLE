#!/usr/bin/env julia
ostr[1] = "   e1=0.00 e2=59.000 e3=6.0 x=5.5678 y=4.56677 state=1 movie=intg1.mrc "
ostr[2] = " e1=5.00 e2=178.000 e3=123.0 x=0.4578 y=2.56677 state=0 movie=intg2.mrc   "






os = Array{typeof(Dict())}()

for i in 1:2
   ostr_split = split(strip(ostr[1]))



end

ostr1_split =

ostr1_dict  = Dict()
for i in eachindex(ostr1_split)
   # println(ostr1_split[i])
   keyval = split(ostr1_split[i], "=")
   ostr1_dict[keyval[1]] = keyval[2]
end
show(ostr1_dict)


# show(ostr_split = split(ostr1, strip(ostr1)))

# ostr2_split = split(strip(ostr1))
