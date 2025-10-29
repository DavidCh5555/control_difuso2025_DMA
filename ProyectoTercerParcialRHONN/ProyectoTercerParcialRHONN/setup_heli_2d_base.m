% SETUP_LAB_HELI_2D
%
% 2 DOF Helicopter (2DHELI) Control Lab: 
% Design of a FF+LQR+I position controller
% 
% SETUP_LAB_HELI_2D sets the model parameters and set the controller
% parameters for the Quanser 2DOF Helicopter system.
%
% Copyright (C) 2005 Quanser Consulting Inc.
% Quanser Consulting Inc.
%
clear all;
%
% ############### USER-DEFINED 2DOF HELI CONFIGURATION ###############
% Cable Gain used for yaw and pitch axes.
K_CABLE_P = 5;
K_CABLE_Y = 3;
% UPM Maximum Output Voltage (V): YAW has UPM-15-03 and PITCH has UPM-24-05
VMAX_UPM_P = 24;
VMAX_UPM_Y = 15;
% Digital-to-Analog Maximum Voltage (V): set to 10 for Q4/Q8 cards
VMAX_DAC = 10;
% Pitch and Yaw Axis Encoder Resolution (rad/count)
K_EC_P = - 2 * pi / ( 4 * 1024 );
K_EC_Y = 2 * pi / ( 8 * 1024 );
% Initial Angle of Pitch (rad)
theta_0 = -40.5*pi/180;
%
% ############### END OF USER-DEFINED DOF HELI CONFIGURATION ###############
%
%
% ############### USER-DEFINED CONTROLLER/FILTER DESIGN ###############
% Anti-windup: integrator saturation (V)
SAT_INT_ERR_PITCH = 5;
SAT_INT_ERR_YAW = 5;
% Anti-windup: integrator reset time (s)
Tr_p = 1;
Tr_y = 1;

% Specifications of a second-order low-pass filter
wcf = 2 * pi * 20; % filter cutting frequency
zetaf = 0.85;        % filter damping ratio
% ############### END OF USER-DEFINED CONTROLLER DESIGN ###############
%
% ############### USER-DEFINED JOYSTICK SETTINGS ###############
% Joystick input X sensitivity used for yaw (deg/s/V)
K_JOYSTICK_X = 85;
% Joystick input Y sensitivity used for pitch (deg/s/V)
K_JOYSTICK_Y = 85;
% Joystick input X sensitivity used for yaw (V/s/V)
K_JOYSTICK_V_X = 10;
% Joystick input Y sensitivity used for pitch (V/s/V)
K_JOYSTICK_V_Y = 5;
% Pitch integrator saturation of joystick (deg)
INT_JOYSTICK_SAT_LOWER = theta_0 * 180 / pi;
INT_JOYSTICK_SAT_UPPER = abs(theta_0)  * 180 / pi;
% Deadzone of joystick: set input ranging from -DZ to +DZ to 0 (V)
JOYSTICK_X_DZ = 0.25;
JOYSTICK_Y_DZ = 0.25;
% ############### END OF USER-DEFINED JOYSTICK SETTINGS ###############
%
% Set the model parameters of the 2DOF HELI.
% These parameters are used for model representation and controller design.
[ K_pp, K_yy, K_yp, K_py, J_eq_p, J_eq_y, B_p, B_y, m_heli, l_cm, g] = setup_heli_2d_configuration();
%
% For the following state vector: X = [ theta; psi; theta_dot; psi_dot]
% Sampling time
Ts=0.001;


