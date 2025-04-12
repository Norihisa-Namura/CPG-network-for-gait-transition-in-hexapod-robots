# Physical simulation code for gait transitions in hexapod robots
# 
# This code requires "hexapod.xml" for simulations
#
# Please run by
#
# mjpython ./mujoco_hexapod.py (for mac)

from scipy.io import loadmat, savemat
import numpy as np
import time
import mujoco

#############
# load data #
#############
gaits = ['wave', 'tetrapod', 'tripod']
name = gaits[0] + '_' + gaits[1] + '_' + gaits[2]
#name = gaits[2] + '_' + gaits[1] + '_' + gaits[0]
mat_data_name = './data/fhn_' + name + '.mat'
mat_data = loadmat(mat_data_name)
ref_angle = mat_data['ref_angle']
l = ref_angle.shape[1]
timestep = mat_data['dt'][0, 0]

################
# mujoco model #
################
model = mujoco.MjModel.from_xml_path('./hexapod.xml')
model.opt.timestep = timestep
data = mujoco.MjData(model)

####################
# initial position #
####################
data.qpos[2] = 0.3 # initial height
# [RF, RM, RH, LF, LM, LH]
init_angle = np.pi * np.array([0,0.25,-0.75,0,0.25,-0.75,0,0.25,-0.75,0,-0.25,0.75,0,-0.25,0.75,0,-0.25,0.75]) # radian
data.qpos[-18:] = init_angle

#########################
# foot tip trajectories #
#########################
link_names = ['RF_tibia', 'RM_tibia', 'RH_tibia', 'LF_tibia', 'LM_tibia', 'LH_tibia']
link_ids = []
for link_name in link_names:
    link_ids.append(model.geom(link_name).id)
height_foot_tip = np.zeros((6, l))

##################
# viewer setting #
##################

count = -1
time_control = 1 # seconds
final = (l - 1) * timestep + time_control # seconds
terminal = False

with mujoco.viewer.launch_passive(model, data) as viewer:
    ##################
    # camera setting #
    ##################
    viewer.cam.type = mujoco.mjtCamera.mjCAMERA_FIXED
    viewer.cam.fixedcamid = model.cam('track').id

    ####################
    # start simulation #
    ####################
    # Close the viewer after the final time
    while viewer.is_running() and data.time < final + timestep/2:
        step_start = time.time()
        
        data.ctrl = init_angle
        if data.time > time_control - timestep/2:
            count += 1
            print(count)

            #################
            # control input #
            #################
            if count == l - 1:
                terminal = True
                break
            else:
                data.ctrl = ref_angle[:, count]
            
            #########################
            # foot tip trajectories #
            #########################
            for k, link_id in enumerate(link_ids):
                pos_link_center = data.geom_xpos[link_id].copy()
                if k < 3:
                    offset_local = np.array([0.2, 0.0, 0.0])
                if k >= 3:
                    offset_local = np.array([-0.2, 0.0, 0.0])
                rotation_matrix = data.geom_xmat[link_id].reshape(3, 3)
                pos_foot_tip = pos_link_center + rotation_matrix @ offset_local
                height_foot_tip[k, count] = pos_foot_tip[2]

        mujoco.mj_step(model, data) # time evolution
        viewer.sync()

        # Time keeping for real-time simulations
        time_until_next_step = model.opt.timestep - (time.time() - step_start)
        if time_until_next_step > 0:
            time.sleep(time_until_next_step)

if terminal:
    save_name = './data/fhn_' + name + '_height_foot_tip.mat'
    dict = {}
    dict['height_foot_tip'] = height_foot_tip

    savemat(save_name, dict)