using PyPlot
using DelimitedFiles
using PyCall

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.08
py = 0.12
sx = 0.83
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(10,5.5))
ax = plt.axes([px,py,sx,sy])
cax = plt.axes([px+sx+dx,py,0.01,sy])

dat = readdlm("GF_DOS_N_16_tso_0.001_tPhU_-0.004_tPhD_-0.004_tsoPh_0.0001_beta_40.0_tl_0.004_V_0.52_P_0.5_coupled.dat")
cp = plt.cm.jet
cmap = plt.cm.ScalarMappable(cmap=cp)
vmn = 0
vmx = 160
ax.imshow(dat[:,2:end],aspect="auto",vmin=vmn,vmax=vmx,interpolation="none",extent=[0,1.250000000000000000e5*fac_t*1e-3,dat[1,1]*1e3,dat[end,1]*1e3],cmap="jet")
ax.set_xlabel(L"$t\,[\mathrm{ps}]$")
ax.set_ylabel(L"$\varepsilon \, [\mathrm{meV}]$",labelpad=-0.5)
cbar = fig.colorbar(cmap, ax=ax, cax = cax, aspect=90,pad=.001)
cbar.set_label(L"$\mathrm{DOS}(\varepsilon)$")
cbar.ax.tick_params(labelsize=fs)


plt.savefig("DOS_time_dep_avg_N_16_tso_0.001_tPhU_-0.004_tPhD_-0.004_tsoPh_0.0001_beta_40.0_tl_0.004_V_0.52_P_0.5.pdf")
