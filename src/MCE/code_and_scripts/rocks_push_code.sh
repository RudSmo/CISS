rm -rf src/__pycache__/
rsync -vax src stefanarne.graeber@file.theorie.physik.uni-goettingen.de:/rocks/stefanarne.graeber/
rsync -vax data_manager.py stefanarne.graeber@file.theorie.physik.uni-goettingen.de:/rocks/stefanarne.graeber/
rsync -vax mce_simulate.py stefanarne.graeber@file.theorie.physik.uni-goettingen.de:/rocks/stefanarne.graeber/
rsync -vax jobscripts stefanarne.graeber@file.theorie.physik.uni-goettingen.de:/rocks/stefanarne.graeber/
#rsync -vax normalize_export_folder_names.py stefanarne.graeber@file.theorie.physik.uni-goettingen.de:/rocks/stefanarne.graeber/
#rsync -vax Makefile stefanarne.graeber@file.theorie.physik.uni-goettingen.de:/rocks/stefanarne.graeber/
