%%%% ECE 4784 Phase 1
% Daniel Russell


%% Time Constants
simTime = 100;% sim run time
hstep = .01; %step size for Eulers method
t = 0:hstep:simTime; 
tau = length(t);

%%
%% Current 
I = zeros(tau); %Current // ONLY USE THIS ONE FOR RESTING MEMBRANE
%I(1:50) = 5; % For a impulse uncomment this. (This gives a
%.5milisecond impulse at 5 microamps/square cm. If you want to increase
%this value change the number after the equals sign.
%I = ones(tau).*10; % this creates a constant current at a given value
%you can change the value at which it is stimulated by changing the 
% number the vector is being mulitplied by.

%% Constants
gK = 36; %max conductance of potassium
gNa = 120; %max conductance of sodium
gL = .3; %max conductance of leakage
EK = -12; %Reversal potential of potassium
ENa = 115; % reversal potential of sodiumm
EL = 10.6; % reversal potential of leakage current
C=1; % cap in microFarads
%Rest potential
Vrest = -70;

%% Gating Variables
% initial states
V = 0;
am = .1*((25-V)/(exp((25-V)/10)-1));
Bm = 4*exp(-V/18);
an = .01*((10-V)/(exp((10-V)/10)-1));
Bn = .125*exp(-V/80);
ah = .07*exp(-V/20);
Bh = (1/(exp((30-V)/10)+1));

% this lets us calculate inital values for m,n, and h 
% initial states

m(1) = am/(am+Bm);
n(1) = an/(an+Bn);
h(1) = ah/(ah+Bh);


%% These are the Currents that are used

%INa = m.^3.*gNa.*h.*(V-ENa);
%IK = n.^4*gK*(V-EK);
%IL = gL*(V-EL);
%Iion = I-IK-INa-IL; %total current

for i=1:tau-1
%% Gated Variables at each time step
am(i) = .1*((25-V(i))/(exp((25-V(i))/10)-1));
Bm(i) = 4*exp(-V(i)/18);
an(i) = .01*((10-V(i))/(exp((10-V(i))/10)-1));
Bn(i) = .125*exp(-V(i)/80);
ah(i) = .07*exp(-V(i)/20);
Bh(i) = (1/(exp((30-V(i))/10)+1));    

%% Currents at each time step
INa = (m(i).^3)*gNa*h(i)*(V(i)-ENa);
IK = (n(i)^4)*gK*(V(i)-EK);
IL = gL*(V(i)-EL);
Iion = I(i)-IK-INa-IL;
%% Calculate derivatives using Euler's    
V(i+1) = V(i) + hstep*Iion/C;
n(i+1) = n(i) + hstep*((an(i) *(1-n(i)) - Bn(i)*n(i)));
m(i+1) = m(i) + hstep*((am(i) *(1-m(i)) - Bm(i)*m(i)));
h(i+1) = h(i) + hstep*((ah(i) * (1-h(i)) -Bh(i)*h(i)));

end
% Puts Plot at the proper starting point
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
title('gK and gNa')
%%







