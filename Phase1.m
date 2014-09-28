%%%% ECE 4784 Phase 1
% Daniel Russell

%% Simulation Parameters

% Constants
simTime = 100;% sim run time
hstep = .01; %step size for Eulers method
t = 0:hstep:simTime;
%I(1:numel(t)) = 0;
tau = length(t);

I = zeros(tau);
gK = 36; %max conductance of potassium
gNa = 120; %max conductance of sodium
gL = .3; %max conductance of leakage
EK = -12; %Reversal potential of potassium
ENa = 115; % reversal potential of sodium
EL = 10.6; % reversal potential of leakage current

C=1; % cap in microFarads

Vrest = -70;

%% Gating Variables
% initial states
V = 0;
am = .1*((25-V)/(exp(25-V/10)-1));
Bm = 4*exp(-V/18);
an = .01*((10-V)/(exp(10-V/10)-1));
Bn = .125*exp(-V/80);
ah = .07*exp(-V/20);
Bh = (1/(exp(30-V/10)+1));

% this lets us calculate inital values for m,n, and h

m(1) = am/(am+Bm);
n(1) = an/(an+Bn);
h(1) = ah/(ah+Bh);


%% Currents

%INa = m.^3.*gNa.*h.*(V-ENa);
%IK = n.^4*gK*(V-EK);
%IL = gL*(V-EL);

%Iion = I-IK-INa-IL; %total current

for i=1:tau-1
%% gated variable at each time steps
am(i) = .1*((25-V(i))/(exp(25-V(i)/10)-1));
Bm(i) = 4*exp(-V(i)/18);
an(i) = .01*((10-V(i))/(exp(10-V(i)/10)-1));
Bn(i) = .125*exp(-V(i)/80);
ah(i) = .07*exp(-V(i)/20);
Bh(i) = (1/(exp(30-V(i)/10)+1));    

I(i) = 0;
%% Currents at each time step
INa = (m(i)^3)*gNa*h(i)*(V(i)-ENa);
IK = (n(i)^4)*gK*(V(i)-EK);
IL = gL*(V(i)-EL);
Iion = I(i) - IK-INa-IL;
%% Calculate derivatives using Euler's    
V(i+1) = V(i) + hstep*Iion/C;
n(i+1) = n(i) + hstep*((an(i) *(1-m(i)) - Bn(i)*n(i)));
m(i+1) = m(i) + hstep*((am(i) *(1-m(i)) - Bm(i)*m(i)));
h(i+1) = h(i) + hstep*((ah(i) * (1-h(i)) -Bh(i)*h(i)));

end
V = V+Vrest;
%% plot
plot(t,V, 'color' ,'r');
hold on
legend({'Voltage'});

ylabel('Voltage (mV)')
xlabel('time (ms)')
title('Voltage vs. Time in Simulated Neuron')

figure
p1 = plot(t,gK.*n.^4,'color','r' , 'linewidth' ,2);
hold on
p2 = plot(t,gNa*(m.^3).*h,'color','g','linewidth' ,2);
%legend([p1,p2], 'Conductance for K', 'Conductance for Na')
ylabel('Conductance')
xlabel('Time (ms)')









