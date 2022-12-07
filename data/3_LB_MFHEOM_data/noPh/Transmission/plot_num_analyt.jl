using LinearAlgebra
using PyPlot
plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)


function TransmissionAnalytical(E1,E2,G,en,P)
           Delta = sqrt.(((en .- E1) .* (en .- E2) .- (G/2 * (1 - (P/2)^2 ))^2 ) .^ 2 .+ (G/2 .* (2 .* en .- (E1 + E2) .+ P/2 * (E1 - E2))) .^ 2)
           Exp1 = (1 + P) .* ( (en .- E2) .^2 .+ (G/2 * ( 1 - P/2))^2 )
           Exp2 = (1 - P) .* ( (en .- E1) .^2 .+ (G/2 * ( 1 + P/2))^2 )
           return G^2 ./ (4 .* Delta .* conj(Delta)) .* (Exp1 .+ Exp2)
end

function TransmissionGreenFunction(E1,E2,G,en,P)
           H = [E1 0;
               0 E2]
           T = []
           GL = G/2 * [1+P 0;
                       0 1-P]
           GR = G/2 * [1 0;
                       0 1]
           Gamma = GL + GR
           for e in en
               Gr = inv(e * Matrix(I,2,2) - H - 1im/2 * Gamma)
               trans = tr(Gr * GL * Gr' * GR)
               push!(T,trans)
           end
           return T
end

E1 = 1.0
E2 = -E1
en = Array(-4:0.01:4.5)
P = 1/2
G = 1/4

fig = plt.figure(figsize = (10,4.5))
px = 0.08
py = 0.12
sx = 0.4
dx = 0.1
sy = 0.85
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px+sx+dx,py,sx,sy])

ax1.plot(en,TransmissionGreenFunction(E1,E2,G,en,P),color="k",linewidth=3,label="Numerical")
ax1.plot(en,TransmissionAnalytical(E1,E2,G,en,P),color="red",linewidth=2,linestyle="dashed",label="Analytical")
ax1.annotate("(a)",(0.05,0.92),xycoords="axes fraction",fontsize=fs)
ax1.set_xlabel(L"$\varepsilon$")
ax1.set_ylabel(L"$T(\varepsilon, P = 0.5)$")
ax1.set_xlim(-4,4.5)
ax1.legend(loc="upper right")


ax2.plot(en, TransmissionAnalytical(E1,E2,G,en,0.5) .- TransmissionAnalytical(E1,E2,G,en,0.0), color="k",linewidth=3)
ax2.annotate("(b)",(0.05,0.92),xycoords="axes fraction",fontsize=fs)
ax2.fill_between(en, TransmissionAnalytical(E1,E2,G,en,0.5) .- TransmissionAnalytical(E1,E2,G,en,0.0), y2 = 0, color = "gray", alpha = 0.6)
#ax2.set_ylim(0.5,1.45)
ax2.set_ylim(-0.2,0.125)
ax2.set_xlabel(L"$\varepsilon$")
ax2.set_ylabel(L"T(\varepsilon, P = 0.5) - T(\varepsilon, P = 0)")
ax2.set_xlim(-4,4.5)

plt.savefig("Comparison_Transmission_Analytical_Numerical.pdf")



