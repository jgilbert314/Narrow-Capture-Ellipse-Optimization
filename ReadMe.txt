ReadMe:

Scripts used in current workflow/testing:
	- TestingDiffTraps
		- Runs optimization for an arrangement of traps
	- CrudeLoop
		- Loops through TestingDiffTraps, changing which traps are the larger ones on each iteration.
		
	- RoughBench
		- A few shoddy scripts used for analyzing failed configurations
	
	- LogSearch
		- Gets index of file with smallest merit function, above a user-defined threshold.
		

Workflow:
	- Run CrudeLoop or TestingDiffTraps
	- Program will run until optimization is complete
	- Suspectedly unnacceptable configurations will be written to file 'TestLog.log' from function 'meritFunctionGeneral'. Conditions specifying log behaviour can be defined here.
	- I usually run CrudeLoop until I see data in the log file (usually takes ~30 sec), then stop the program and analyze.
	- Run LogSearch to find optimum arrangement and trap indices. Plot configuration using PlotBenchEllipse.
	- Currently the program runs without catching bad configurations, but inspecting the results shows that the results are unnacceptable.