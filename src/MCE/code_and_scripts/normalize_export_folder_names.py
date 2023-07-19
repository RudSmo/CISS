import os
import sys

if len(sys.argv) == 1:
	folders_to_check = [ f.path for f in os.scandir("data_export/") if f.is_dir() ]
else:
	#folders_to_check = [ f.replace("\\", "").replace("(", r"\(").replace(")", r"\)") for f in sys.argv[1:]]
	folders_to_check = [ f for f in sys.argv[1:] if os.path.isdir(f)]

print_changes = True
apply_changes = True

temp_folder = "data_export/temp/"	# must have a trailing "/"

if os.path.isdir(temp_folder) and len([ f.path for f in os.scandir("data_export/")]) != 0:
	print(f"{temp_folder} must not exist or be empty when starting this script!")
	quit(1)
os.system("mkdir -p " + temp_folder)


for d in folders_to_check:
	
	if print_changes:
		print(f"\n{d}:")
	subfolders = [f.path for f in os.scandir(d) if f.is_dir() ]
	subfolders = [f[len(d):] for f in subfolders]
	subfolders.sort()
	
	wrong_placed_subfolders = []
	for i in range(len(subfolders)):
		sf = d + subfolders[i]
		for f in os.scandir(sf):
			if f.is_dir():
				wrong_placed_subfolders.append(f.path[len(d):])
		
	#print(f"{wrong_placed_subfolders=}")
	subfolders = wrong_placed_subfolders + subfolders
	
	d = d.replace("(", r"\(").replace(")", r"\)")
	for i in range(len(subfolders)):
		sf = d + subfolders[i]		
		sf_new = temp_folder + str(i)
		mv_command = f"mv {sf} {sf_new}"
		if print_changes:
			print(mv_command)
		if apply_changes:
			os.system(mv_command)
	
	mv_command = f"mv -t {d} {temp_folder}*"
	if print_changes:
		print(mv_command)
	if apply_changes:
		os.system(mv_command)

os.system("rmdir " + temp_folder)
