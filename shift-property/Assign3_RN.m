%% SECTION 1 - Plotting the Shift Property (Logarathmic Increase)  
%For four scenarios, four lists to store the activity of the neurons
x1 = [];
x2 = [];
x3 = [];
x4 = [];
%Decay term
A = 1;
%Upper bound
B = 10;
%Total time for which input will be increasing
T = 1000000;
%Timestep
dt = 0.1;
%Inputs will increase by dt until they reach T
I = 0:dt:T;
%Set M to be log(I) 
M = log(I);
%Inhibitory inputs from the other neurons for each of the four scenarios
L1 = 10;
L2 = 100;
L3 = 1000;
L4 = 10000;
%set y to be the length of the input array to iterate through more easily
y = max(size(I));


%iterate through the number of inputs to find the neuron activity in each 
%of the four scenarios. this equation is based on solving at equilibrium
for i=1:y
   x1(i) = (B * exp(M(i))) / (A + exp(M(i)) + L1);
   x2(i) = (B * exp(M(i))) / (A + exp(M(i)) + L2);
   x3(i) = (B * exp(M(i))) / (A + exp(M(i)) + L3);
   x4(i) = (B * exp(M(i))) / (A + exp(M(i)) + L4);
end

%figure 1
figure;
plot(M,x1,'LineWidth',3)
xlabel("log(I)")
ylabel("Activity")
title(["The Sensitivity of Neurons with Logarithmically Varying Background ",
    "Intensities in a Feedforward Shunting Network"]);

hold on
plot(M,x2,'LineWidth',3)

hold on
plot(M,x3,'LineWidth',3)

hold on
plot(M,x4,'LineWidth',3)

legend("L=10","L=100","L=1000","L=10000")


%% SECTION 2 - Plotting the Shift Property (Linear Increase)  

%Inhibitory inputs from the other neurons for each of the four scenarios
L1 = 10;
L2 = 20;
L3 = 30;
L4 = 40;

%iterate through the number of inputs to find the neuron activity in each 
%of the four scenarios. this equation is based on solving at equilibrium
for i=1:y
   x1(i) = (B * exp(M(i))) / (A + exp(M(i)) + L1);
   x2(i) = (B * exp(M(i))) / (A + exp(M(i)) + L2);
   x3(i) = (B * exp(M(i))) / (A + exp(M(i)) + L3);
   x4(i) = (B * exp(M(i))) / (A + exp(M(i)) + L4);
end

%figure 2
figure;
plot(M,x1,'LineWidth',3)
xlabel("log(I)")
ylabel("Activity")
title(["The Sensitivity of Neurons with Linearly Varying Background ",
    "Intensities in a Feedforward Shunting Network"]);

hold on
plot(M,x2,'LineWidth',3)

hold on
plot(M,x3,'LineWidth',3)

hold on
plot(M,x4,'LineWidth',3)

legend("L=10","L=20","L=30","L=40")

%% SECTION 3 - Plotting Time vs Neuron Activity

figure;

%create and plot array of inputs
a = zeros(1,100);
b = ones(1,100) * 1;
I = [a b a b b a b b b a];
plot(I,'LineWidth',2)

%total run time
time = 10;
%timestep
dt = 0.01;
%initalize array of activity values, which is the y-axis of plots
x = [0];
%array of the derivative at each t value
dxdt = [];

L = .00001;
%get x values for each time through Euler's method
for t=0+dt:dt:time
    %position = index of the array, convert to integer
    pos = int16(t/dt);
    %find derivative, multiply by timestep, add to last known value
    dxdt(pos) = (-A * x(pos)) + ((B - x(pos)) * I(pos)) - (x(pos) * L);
    x(pos+1) = x(pos) + (dt*dxdt(pos));
end 
hold on
plot(x,'LineWidth',2)

L = 1;
for t=0+dt:dt:time
    %position = index of the array, convert to integer
    pos = int16(t/dt);
    %find derivative, multiply by timestep, add to last known value
    dxdt(pos) = (-A * x(pos)) + ((B - x(pos)) * I(pos)) - (x(pos) * L);
    x(pos+1) = x(pos) + (dt*dxdt(pos));
end 
hold on
plot(x,'LineWidth',2)

L = 10;
%get x values for each time through Euler's method
for t=0+dt:dt:time
    %position = index of the array, convert to integer
    pos = int16(t/dt);
    %find derivative, multiply by timestep, add to last known value
    dxdt(pos) = (-A * x(pos)) + ((B - x(pos)) * I(pos)) - (x(pos) * L);
    x(pos+1) = x(pos) + (dt*dxdt(pos));
end 
hold on
plot(x,'LineWidth',2)

legend("Input", "L=0.1", "L=1","L=10");
xlabel("Time")
ylabel("Activity")
title(["Activity of a Neruon in a Feedforward Shunting Network", ...
    "Over Time with Varying Levels of Background"])

%% SECTION 4 - Plotting Time vs Neuron Activity Plus C Term 
figure();
a = zeros(1,100);
b = ones(1,100) * 1;
I = [a b a b b a b b b a];
plot(I,'LineWidth',2)

%total run time
time = 10;
%timestep
dt = 0.01;
%array of time values which is the x-axis of plots
times = 0:time/dt;
%begin array of activity values, which is the y-axis of plots, with init
x = [0];
%array of the derivative at each t value
dxdt = [];

%%Part a
A = 1;
B = 10;
C = 2;
L = .1;

%get x values for each time through Euler's method
for t=0+dt:dt:time
    %position = index of the array, convert to integer
    pos = int16(t/dt);
    %find derivative, multiply by timestep, add to last known value
    dxdt(pos) = (-A * x(pos)) + ((B - x(pos)) * I(pos)) - ((x(pos) + C) * L);
    x(pos+1) = x(pos) + (dt*dxdt(pos));
end 
hold on
plot(x,'LineWidth',2)

L = 1;
for t=0+dt:dt:time
    %position = index of the array, convert to integer
    pos = int16(t/dt);
    %find derivative, multiply by timestep, add to last known value
    dxdt(pos) = (-A * x(pos)) + ((B - x(pos)) * I(pos)) - ((x(pos) + C) * L);
    x(pos+1) = x(pos) + (dt*dxdt(pos));
end 
hold on
plot(x,'LineWidth',2)

L = 10;

for t=0+dt:dt:time
    %position = index of the array, convert to integer
    pos = int16(t/dt);
    %find derivative, multiply by timestep, add to last known value
    dxdt(pos) = (-A * x(pos)) + ((B - x(pos)) * I(pos)) - ((x(pos) + C) * L);
    x(pos+1) = x(pos) + (dt*dxdt(pos));
end 
hold on
plot(x,'LineWidth',2)


xlabel("Time")
ylabel("Activity")
title(["Activity of a Neruon in a Feedforward Shunting Network", ...
    "Over Time with Varying Levels of Background and a Hyperpolarization Term"])
legend("Input", "L=0.1", "L=1","L=10");

%% SECTION 5 - Plotting Varying C Terms
figure();
a = zeros(1,100);
b = ones(1,100) * 1;
I = [a b a b b a b b b a];
plot(I,'LineWidth',2)

%total run time
time = 10;
%timestep
dt = 0.01;
%array of time values which is the x-axis of plots
times = 0:time/dt;
%begin array of activity values, which is the y-axis of plots, with init
x = [0];
%array of the derivative at each t value
dxdt = [];

%%Part a
A = 1;
B = 10;
C = 0.01;
L = 0.1;

%get x values for each time through Euler's method
for t=0+dt:dt:time
    %position = index of the array, convert to integer
    pos = int16(t/dt);
    %find derivative, multiply by timestep, add to last known value
    dxdt(pos) = (-A * x(pos)) + ((B - x(pos)) * I(pos)) - ((x(pos) + C) * L);
    x(pos+1) = x(pos) + (dt*dxdt(pos));
end 
hold on
plot(x,'LineWidth',2)

C = 1;
for t=0+dt:dt:time
    %position = index of the array, convert to integer
    pos = int16(t/dt);
    %find derivative, multiply by timestep, add to last known value
    dxdt(pos) = (-A * x(pos)) + ((B - x(pos)) * I(pos)) - ((x(pos) + C) * L);
    x(pos+1) = x(pos) + (dt*dxdt(pos));
end 
hold on
plot(x,'LineWidth',2)

C = 10;
for t=0+dt:dt:time
    %position = index of the array, convert to integer
    pos = int16(t/dt);
    %find derivative, multiply by timestep, add to last known value
    dxdt(pos) = (-A * x(pos)) + ((B - x(pos)) * I(pos)) - ((x(pos) + C) * L);
    x(pos+1) = x(pos) + (dt*dxdt(pos));
end 
hold on
plot(x,'LineWidth',2)


xlabel("Time")
ylabel("Activity")
title(["Activity of a Neruon in a Feedforward Shunting Network", ...
    "Over Time with Varying Hyperpolarization Terms"])
legend("Input", "C=0.01", "C=1","C=10");
%% SECTION 6 - Plotting Time vs Neuron Activity Plus C Term in terms of I vs A
figure();
a = zeros(1,100);
b = ones(1,100) * 1;
I = [a b a b b a b b b a];
plot(I,'LineWidth',2)

%total run time
time = 10;
%timestep
dt = 0.01;
%array of time values which is the x-axis of plots
times = 0:time/dt;
%begin array of activity values, which is the y-axis of plots, with init
x = [0];
%array of the derivative at each t value
dxdt = [];

% I > A
A = .01;
B = 10;
C = 1;
L = .01;

%get x values for each time through Euler's method
for t=0+dt:dt:time
    %position = index of the array, convert to integer
    pos = int16(t/dt);
    %find derivative, multiply by timestep, add to last known value
    dxdt(pos) = (-A * x(pos)) + ((B - x(pos)) * I(pos)) - ((x(pos) + C) * L);
    x(pos+1) = x(pos) + (dt*dxdt(pos));
end 
hold on
plot(x,'LineWidth',2)

% I < A
L = 2;
A = 5;
for t=0+dt:dt:time
    %position = index of the array, convert to integer
    pos = int16(t/dt);
    %find derivative, multiply by timestep, add to last known value
    dxdt(pos) = (-A * x(pos)) + ((B - x(pos)) * I(pos)) - (x(pos) * C * L);
    x(pos+1) = x(pos) + (dt*dxdt(pos));
end 
hold on
plot(x,'LineWidth',2)

%I ~= A
L = 1;
A = 2;
for t=0+dt:dt:time
    %position = index of the array, convert to integer
    pos = int16(t/dt);
    %find derivative, multiply by timestep, add to last known value
    dxdt(pos) = (-A * x(pos)) + ((B - x(pos)) * I(pos)) - (x(pos) * C * L);
    x(pos+1) = x(pos) + (dt*dxdt(pos));
end 
hold on
plot(x,'LineWidth',2)


xlabel("Time")
ylabel("Activity")
title(["Activity of a Neruon in a Feedforward Shunting Network", ...
    "Over Time with Varying Relationships Between I and A"])
legend("Input", "I > A", "I < A","I~=A","Noise");

%% SECTION 7 - Plotting the Shift Property with a Small Hyperpolarization Term
A = 1;
B = 10;
C = .1;

%Inhibitory inputs from the other neurons for each of the four scenarios
L1 = 10;
L2 = 100;
L3 = 1000;
L4 = 10000;

%iterate through the number of inputs to find the neuron activity in each 
%of the four scenarios. this equation is based on solving at equilibrium
for i=1:y
   x1(i) = ((B * exp(M(i))) - (C * L1)) / (A + exp(M(i)) + L1);
   x2(i) = ((B * exp(M(i))) - (C * L2)) / (A + exp(M(i)) + L2);
   x3(i) = ((B * exp(M(i))) - (C * L3)) / (A + exp(M(i)) + L3);
   x4(i) = ((B * exp(M(i))) - (C * L4)) / (A + exp(M(i)) + L4);
end

%figure 1
figure;
plot(M,x1,'LineWidth',3)
xlabel("log(I)")
ylabel("Activity")
title(["The Sensitivity of Neurons with Linearly Varying Background Intensities",
    "in a Feedforward Shunting Network with a Low Hyperpolarization Term"]);

hold on
plot(M,x2,'LineWidth',3)

hold on
plot(M,x3,'LineWidth',3)

hold on
plot(M,x4,'LineWidth',3)

legend("L=10","L=100","L=1000","L=10000")

%% SECTION 8 - Plotting the Shift Property with a Large Hyperpolarization Term
C = 100;

%Inhibitory inputs from the other neurons for each of the four scenarios
L1 = 10;
L2 = 100;
L3 = 1000;
L4 = 10000;

%iterate through the number of inputs to find the neuron activity in each 
%of the four scenarios. this equation is based on solving at equilibrium
for i=1:y
   x1(i) = ((B * exp(M(i))) - (C * L1)) / (A + exp(M(i)) + L1);
   x2(i) = ((B * exp(M(i))) - (C * L2)) / (A + exp(M(i)) + L2);
   x3(i) = ((B * exp(M(i))) - (C * L3)) / (A + exp(M(i)) + L3);
   x4(i) = ((B * exp(M(i))) - (C * L4)) / (A + exp(M(i)) + L4);
end

%figure 1
figure;
plot(M,x1,'LineWidth',3)
xlabel("log(I)")
ylabel("Activity")
title(["The Sensitivity of Neurons with Linearly Varying Background Intensities",
    "in a Feedforward Shunting Network with a Large Hyperpolarization Term"]);

hold on
plot(M,x2,'LineWidth',3)

hold on
plot(M,x3,'LineWidth',3)

hold on
plot(M,x4,'LineWidth',3)

legend("L=10","L=100","L=1000","L=10000")

