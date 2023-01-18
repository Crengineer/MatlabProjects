%% Controller Parameters

kd = 100;
ki = 50;
kp = 500;

%% Quarter-Car model Parameters

m = 1;
b = 10;
k = 20;

%% Simulation Parameters

Sim_time= 7;
Sim_step = 1;

%% Simulink execution

sim('Project8Simulink');
time = out.x.time;
pos  = out.x.data;
step_resp = out.step.data;

figure(1)
plot(time,pos);
hold all;
grid on;
xlabel('Time');
ylabel('Position');
plot(time,step_resp);

