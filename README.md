# <b> A central pattern generator network for simple control of gait transitions in hexapod robots based on phase reduction </b>

## <b> \<Description\> </b>
<p>
The code is for simulating the gait transitions in hexapod robots by the central pattern generator network proposed in the paper:
</p>

<p>
N. Namura and H. Nakao, 
<a href="https://link.springer.com/article/10.1007/s11071-024-10773-x" target="_blank">
"A central pattern generator network for simple control of gait transitions in hexapod robots based on phase reduction,"
</a>
<i> Nonlinear Dynamics </i> <b> 113</b>, 10105–10125 (2025).
</p>

<p>
If you use this code, please cite this paper.
</p>


## <b> \<Usage\> </b>
For simulating gait transitions, please run the codes by the following procedure:
1. Run "main.m" in MATLAB.
2. Run "mujoco_hexapod.py" by python ("mjpython" should be used instead of "python" in Mac environment for running this code).
3. If "show_results.m" is run after the python code, you can find orange boxes in panel (b).


### <b> Files: </b>
- fitzhugh_nagumo.m: FitzHugh–Nagumo (FHN) oscillator
- funcs.m: Set of functions
- hexapod.xml: XML file for modeling the hexapod robot in physical simulation
- main.m: Main script for the simulations in MATLAB
- mujoco_hexapod.py: Script for physical simulation and visualization of the gait transitions by 
<a href="https://github.com/google-deepmind/mujoco" target="_blank">
MuJoCo
</a>
- mystyle.m: Settings for MATLAB
- PSF.m: Calculation of the PSF
- reference_trajectory.m: Calculation of reference joint angles from CPG outputs
- show_results.m: Visualization for Fig. 8 and 10
- utils.m: Utility functions mainly for visualizing


## <b> \<Environments\> </b>
### OS:
- MacOS Sonoma 14.5
  
### MATLAB R2022a
- No specific toolbox is required

### Python 3.9.19
- requirements.txt is attached


## <b> \<Contact\> </b>
E-mail: namura.n.aa@m.titech.ac.jp