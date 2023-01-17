%% Data taken from .xlsx file and vectors created

table = readtable('Battery_Parameters.xlsx');
StateOfCharge___ = table(:,1);
OpenCircuitVoltage_V_ = table(:,2);
ChargeResistance_Ohms_ = table(:,3);
DischargeResistance_Ohms_ = table(:,4);
s = max(size(StateOfCharge___));
base_vect = linspace(0,1,s);
%% Figures

figure(1);
grid on;
plot(base_vect',table.StateOfCharge___);

figure(2);
grid on;
plot(base_vect',table.OpenCircuitVoltage_V_);

figure(3);
grid on;
plot(base_vect',table.ChargeResistance_Ohms_);

figure(4);
grid on;
plot(base_vect',table.DischargeResistance_Ohms_);

%% As shown by the course

Data = xlsread('Battery_Parameters.xlsx');
SOC = Data(:,1);
OCV = Data(:,2);
R_charge = Data(:,3);
R_discharge = Data(:,4);
I = 2.3; % Available current
Cn = 2.3 * 3600; %Capacity
Sim_time = 3600; % Simulation time

%% Plots

figure(5);
plot(SOC,OCV);
grid on;
xlabel('SOC');
ylabel('OCV');

figure(6);
plot(SOC,R_charge);
grid on;
xlabel('SOC');
ylabel('R__charge');

figure(7);
plot(SOC,R_discharge);
grid on;
xlabel('SOC');
ylabel('R__discharge');

%% Simulink run

sim('Project6Simulink');

