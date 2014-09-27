%% ECE 4784 Phase 1
% Daniel Russell

%% Simulation Parameters

% Constants

gK = 36;
gNa = 120;
gL = .3;
EK = -12;
ENa = 115;
EL = 10.6;

Vrest = -70;

vm = Vrest;
%% Gating Variables

am = .1*(25-vm)/(e(25-vm/10)-1);
Bm = 4*e(-vm/18);
an = .1*(10-vm)/(e(10-vm/10)-1);
Bn = .125*e(-vm/80);
ah = .7*e(-vm/20);
Bh = .1*(1/(e(30-vm/10)+1));

%% Currents

INa = m^3*gNa*h*(vm-ENa);
IK = n^4*gK*(vm-EK);
IL = gL*(vm-EL);
Iion = I-IK-INa-IL;

%% Derivatives

plot(t,V);
hold on
legend({'Voltage'});

ylabel('Voltage (mV)')
xlabel('time (ms)')
title('Voltage vs. Time in Simulated Neuron')

figure
p1 = plot(t,gK*n.^4);
hold on
p2 = plot(t,gNa*(m.^3).*h);
legend([p1,p2], 'Conductance for K', 'Conductance for Na')
ylabel('Conductance')
xlabel('time (ms)')







