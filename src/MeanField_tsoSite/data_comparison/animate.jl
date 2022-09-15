using Plots
using DelimitedFiles


dat = readdlm("L40_Hop+Onsite+SOC_dt0.001_1_electrons.dat")
anim = @animate for i = 2:2:length(dat[2:end,1])
    plot(dat[i,2:2:end],legend=false,ylims=(0,1))
end


gif(anim, "N40_fps15.gif", fps = 30)
