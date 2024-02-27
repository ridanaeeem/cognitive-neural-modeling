%input strength
input = 5;
%start time of input
t1 = 1;
%end time of input
t2 = 6;
%total run time
time = 10;
%array of time values which is the x-axis of plots
times = 0:10;
%array of position values which is the y-axis of plots
%it is time+1 to account for matlab indexing at 1
x = zeros(1,time+1);
%array of the derivative at each t value
dxdt = zeros(1,time + 1);

%%Simulation 1
%Decay term
A = 1;

%get x values for each time through Euler's method
for t=1:time
    if t >= t1 && t <= t2
        I = input;
    else
        I = 0;
    end
    
    %x(t) = (-5*exp(-A * t)) + (I/A);
    dxdt(t) = (-A * x(t)) + I;
    x(t+1) = x(t) + dxdt(t);
end

%plot
plot(times,x)
ylim([0,6])

%%Simulation 2
%change the decay term
A = 2;

%get x values for each time through Euler's method
for t=1:time
    if t >= t1+1 && t <= t2+1
        I = input;
    else
        I = 0;
    end
    
    %x(t) = (-5*exp(-A * t)) + (I/A);
    dxdt(t) = (-A * x(t)) + I;
    x(t+1) = x(t) + dxdt(t);
end

%plot
figure;
plot(times,x)
ylim([0,6])



% init = 0;
% x(1) = init;


