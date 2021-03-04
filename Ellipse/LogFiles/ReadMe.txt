This folder contains the files used to find optimal configurations of up to 50 traps, and records of all intermediate configurations. Some intermediate configurations may be invalid. If so it should be obvious, as the merit will be either 0, NaN, or an absurdly large number. 

Folder names: 
Use the format
	Log_eps0p***-ecc0p***
where 'eps' is short for 'epsilon', the size of the traps, and 'ecc' is short for 'eccentricity', the eccentricity of the domain. In each case, '0p' stands for '0.'.

Folder contents:
Each folder should contain the following files
	- WorkbenchEllipse_eps0p***-ecc0p***
		Matlab script used to specify optimization problem and create log files
	- log*.csv
		Running log files. Each iteration of optimization will be saved in a numbered log file.
	- GoodLog***.csv
		Composite log files. After each iteration of optimization the newly computed optimums are compared to the old, and a new composite list is produced which contains the best known results.
The format of the log files is
	N, coordinates, merit function value, computation time
where each field is delimited by a ',' character, and there are 2N coordinates in the order [x1, ..., xN, y1, ..., yN]

