# shell scripts
Here lie bits of shell and python code which are useful for saving, renaming,
and storing data generated during this project. Below are short descriptions of
each function and how they should be used. **For all python scripts,
information of usage can be found by typing `./script.py --help`.**

## `fcs_rename.py`
**language** : python

**purpose** : Easily rename a large group of files generated using MACSQuant
	flow cytometer. 

**arguments** :

	-d, --dir : DIRECTORY
		Name of directory containing files to be renamed.

	-t, --template : CSV FILE
		Path to csv containing new names of files. The structure of
		this csv file should be the same as the layout of samples in
		the 96-well plates. Please stick to the standardized naming
		convention of `YYYYMMDD_strain_operator_rbs_xuMIPTG``

	-e, --ext : PATTERN
		Extension to append to each file. If not specified, '.fcs' will
		be used.

	-o, --output : DIRECTORY
		Path to output directory. If none is provided, files will be
		renamed in place. If this directory does not exist, it will be
		automatically made. If it exists but is not empty, the user
		will be prompted for confirmation that the renaming should
		proceed.

	-v, --verbose : NONE
		Optional flag to print progress of renaming to screen.

	-f, --force : NONE
		Optional flag to force renaming of files is ouput directory is
		not empty.

**example** :

This example assumes that the user is standing in the `/data/flow/fcs/` folder,
and that the `YYYYMMDD_fx_fcs_rename.csv` file that lists how to rename the
files already exists in the same folder. Then to rename the files one just has
to run:

	>$ ../../../code/shell_scripts/fcs_rename.py -d ./ -t 20191022_r1_fcs_rename.csv -p RP2019-10-22 -vf 
		
## `processing_fcs.py`
**language** : python

**purpose** : Read a provided Flow Cytometry Standard (fcs) file or directory
	of fcs files, extract the desired channels, and save them as a Comma Separated
	Value (csv) file in a specified output folder.

**arguments** : 

	-i, --input_file :  FILENAME
		Name of individual file to process.

	-d, --directory : DIRECTORY
		Name of directory containing files to be processed.	

	-p, --pattern : PATTERN
		Pattern within desired filenames to be processed.
	
	-o, --output : DIRECTORY
		Path to output directory. If this directory does not exist, it
		will be made. If this directory exists but is not empty, user
		input will be required to continue with execution.

	-c, --channel : CHANNEL NAME
		Name of channel (exact) to be extracted from the fcs file. If
		multiple channels are desired, multiple calls of -c must be
		called.

	-v, --verbose : NONE
		Optional flag to print progress of processing to screen.

	-f, --force : NONE
		Optional flag to force renaming of files if output directory is
		not empty.

**example** : 

This example assumes that the user is standing in the `/data/flow/fcs/` folder,
and that the `csv` files will be exported to `/data/flow/csv/`. Also the
example assumes that multiple files will be processed, all of them having a
common pattern `20191022_r1_` in the filename. Then all the user needs to write
is:

    >$ ../../../code/shell_scripts/fcs_processing.py -d ./ -p 20191022_r1_ -o ../csv -vf
