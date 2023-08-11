using PyPlot
using LinearAlgebra
using3D()

off = 0.5
theta = Array(-2*pi-off:(4*pi+2*off)/150:2*pi+off)
r = 1.2
z = Array(-2:4/150:2)
x = r*sin.(theta)
y = r*cos.(theta)

num_sites = 3
th_sites = 2*pi/num_sites
fitting = (4*pi)/th_sites
th_coord = [theta[1]+off+i*th_sites for i in 0:fitting]

inds = []
for (i,th) in enumerate(theta)
	for t in th_coord
		if abs(th - t) < 0.45*1e-1
			push!(inds,i)
			break
		end
	end
end

		
xs = r*sin.(th_coord)
ys = r*cos.(th_coord)
zs = z[inds]


xp = [-r-off,r+off]
yp = xp'
zp = [-1,1]#[(zs[num_sites]-zs[num_sites-1])/2, (zs[2*num_sites]-zs[2*num_sites-1])/2]


size = 220
ax = plt.axes(projection="3d")
ax.view_init(elev=14., azim=40)
ax.plot3D(x,y,z,color="k")
ax.scatter3D(xs[num_sites:2*num_sites-1],ys[num_sites:2*num_sites-1],zs[num_sites:2*num_sites-1],color="blue",s=size)
ax.scatter3D(xs[1:num_sites-1],ys[1:num_sites-1],zs[1:num_sites-1],color="gray",s=size)
ax.scatter3D(xs[2*num_sites:end],ys[2*num_sites:end],zs[2*num_sites:end],color="gray",s=size)

ax.set_axis_off()

surf(xp,yp',0 * xp .+ 0 * yp .+ zp[1],color="gray",alpha=0.5,linewidth=0,edgecolors="gray")#,antialiased=false, rcount=200, ccount=200)
surf(xp,yp',0 * xp .+ 0 * yp .+ zp[2],color="gray",alpha=0.5,linewidth=0,edgecolors="gray")#,antialiased=false, rcount=200, ccount=200)
plt.savefig("Helix_3site.pdf")

#plt.show()




