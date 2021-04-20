# Optimal-Regulation-of-Cart-Inverted-Pendulum-System-Comparison-of-HJB-and-Classical-Methodologies

__README__

__Overview of Programs.__ The MATLAB and Simulink programs developed for this project are organized
as follows:

* cip main.m: Main MATLAB program, controls calls to all other programs and Simulink models.
* cip sim.slx: Main Simulink model, controls simulation functions for all designs and test scenarios.
* util: Folder containing all supporting functions (i.e., GMS optimizer, figure saving functions, etc.)

__How to Run the Code.__

1. Open cip main.m.
2. Open cip sim.slx.
3. Select the desired test to conduct via the variable sim type (see line 146). This variable determines the
type of test conducted by the simulator cip sim.slx, and has the following options:

	* Value = 1: Initial condition response.
	* Value = 2: Disturbance force response.
4. Select the desired plant type via the variable p select (see line 219). This variable determines the type
of plant used in the simulator cip sim.slx, and has the following options:
	* Value = 1: Plant linearized about a standing pendulum position.
	* Value = 2: Full nonlinear plant.
5. Click the “Run” button on cip main.m, which will automatically initialize, simulate, and plot the results
for all three design methodologies.
6. Some additional variables to control program settings are listed below:
	* __savefigs__ (line 73): Set to 1 to automatically save all plotted figures. Set to 0 to not save any
figures. Figures will be saved to the folder named figures located in the same directory location
as cip main.m.
	* __plot_individual__ (line 222): Set to 1 to plot simulation results for each of the three methodologies
individually, as well as plot simulation results for all three mothodologies simultaneously on the
same plot. Set to 0 to only plot simulation results for all three mothodologies simultaneously on
the same plot.
	* __plot_freq_resp__ (line 225): Set to 1 to plot frequency response data for the GMS design. Set to 0
to not plot this data.
