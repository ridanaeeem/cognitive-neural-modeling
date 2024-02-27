close all

tStart = 0;
tEnd = 1;
time = 10;
dt = 0.1;
times = [0:time/dt];
times = times * dt;

A = 1;
B = 3;
x1 = zeros(time/dt,10);
x2 = zeros(time/dt,10);
x3 = zeros(time/dt,10);
x4 = zeros(time/dt,10);
I = [ 0.2, 0.6, 0.9, 0.6, 0.2, 0.1, 0.4, 0.8, 0.4, 0.1 ];

for t=0+dt:dt:time
    % input is only presented from t = 0 to 1
    if t >= tStart && t <= tEnd
        %first input
        %I = [ 0.2, 0.6, 0.9, 0.6, 0.2, 0.1, 0.4, 0.8, 0.4, 0.1 ];
        %second input
        I = [0.7, 0.6, 0.8, 0.9, 0.5, 0.3, 0.5, 0.7, 0.8, 0.4];
    else
        I = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
    end

    % go through each neuron in the network
    for i=1:10
        %convert time to an integer
        pos = int16(t/dt);
        F = 0.25;
        %linear
        f1 = x1(pos,i);
        signalExceptI1 = 0;
        for k=1:10
            if k ~= i
                signalExceptI1 = signalExceptI1 + x1(pos,k);
            end
        end
        %faster than linear
        f2 = x2(pos,i) * x2(pos,i);
        signalExceptI2 = 0;
        for k=1:10
            if k ~= i
                signalExceptI2 = signalExceptI2 + (x2(pos,k) * x2(pos,k));
            end
        end
        %slower than linear
        f3 = x3(pos,i) / (F + x3(pos,i));
        signalExceptI3 = 0;
        for k=1:10
            if k ~= i
                signalExceptI3 = signalExceptI3 + ((x3(pos,k) / (F + x3(pos,k))));
            end
        end
        %sigmoid
        f4 = (x4(pos,i) * x4(pos,i)) / (F + (x4(pos,i) * x4(pos,i)));
        signalExceptI4 = 0;
        for k=1:10
            if k ~= i
                signalExceptI4 = signalExceptI4 + (x4(pos,k) * x4(pos,k)) / (F + (x4(pos,k) * x4(pos,k)));
            end
        end
        % change in activity
        dxdt1 = (-A * x1(pos,i)) + ((B - x1(pos,i)) * (f1 + I(i))) - (x1(pos,i) * signalExceptI1);
        dxdt2 = (-A * x2(pos,i)) + ((B - x2(pos,i)) * (f2 + I(i))) - (x2(pos,i) * signalExceptI2);
        dxdt3 = (-A * x3(pos,i)) + ((B - x3(pos,i)) * (f3 + I(i))) - (x3(pos,i) * signalExceptI3);
        dxdt4 = (-A * x4(pos,i)) + ((B - x4(pos,i)) * (f4 + I(i))) - (x4(pos,i) * signalExceptI4);
        % activity
        % row = pos = time (out of time/dt)
        % col = i neuron index 
        y = i;
        x1(pos + 1,i) = x1(pos,i) + (dt * dxdt1);
        x2(pos + 1,i) = x2(pos,i) + (dt * dxdt2);
        x3(pos + 1,i) = x3(pos,i) + (dt * dxdt3);
        x4(pos + 1,i) = x4(pos,i) + (dt * dxdt4);
    end
end
% normalized activity
X1 = x1/sum(sum(x1));
X2 = x2/sum(sum(x2));
X3 = x3/sum(sum(x3));
X4 = x4/sum(sum(x4));

X = [1:10];
Y = times;
%normalized f1 vs time - each with 10 curves
figure;
plot(times,X1,'LineWidth',2)
title("Normalized Activity of Neurons in a Recurrent Neural Network with Signal" + ...
    " Function 1 Over Time")
legend("Neuron 1", "Neuron 2", "Neuron 3",  "Neuron 4", ...
     "Neuron 5",  "Neuron 6",  "Neuron 7",  "Neuron 8", ...
      "Neuron 9",  "Neuron 10");
xlabel("Time")
ylabel("Neuron activity")

figure;
surf(X,Y,X1)
title("Normalized Activity of Neurons in a Recurrent Neural Network with Signal" + ...
    " Function 1 Over Time - Surface Plot")
xlabel("Neuron Index")
ylabel("Time")
zlabel("Neuron activity")

%final normalized f1 vs index
figure;
plot(X1(101, :), "LineWidth",2) 
title("The Final Normalized Activity for Neurons in a Recurrent Neural Network with Signal Function 1")
xlabel("Neuron Index")
ylabel("Neuron Activity")
xlim([1 10])

%normalized f2 vs time - each with 10 curves
figure;
plot(times,X2,'LineWidth',2)
title("Normalized Activity of Neurons in a Recurrent Neural Network with Signal" + ...
    " Function 2 Over Time")
legend("Neuron 1", "Neuron 2", "Neuron 3",  "Neuron 4", ...
     "Neuron 5",  "Neuron 6",  "Neuron 7",  "Neuron 8", ...
      "Neuron 9",  "Neuron 10");
xlabel("Time")
ylabel("Neuron activity")

figure;
surf(X,Y,X2)
title("Normalized Activity of Neurons in a Recurrent Neural Network with Signal" + ...
    " Function 2 Over Time - Surface Plot")
xlabel("Neuron Index")
ylabel("Time")
zlabel("Neuron activity")

%final normalized f2 vs index
figure;
plot(X2(101, :), "LineWidth",2) 
title("The Final Normalized Activity for Neurons in a Recurrent Neural Network with Signal Function 2")
xlabel("Neuron Index")
ylabel("Neuron Activity")
xlim([1 10])

%normalized f3 vs time - each with 10 curves
figure;
plot(times,X3,'LineWidth',2)
title("Normalized Activity of Neurons in a Recurrent Neural Network with Signal" + ...
    " Function 3 Over Time")
legend("Neuron 1", "Neuron 2", "Neuron 3",  "Neuron 4", ...
     "Neuron 5",  "Neuron 6",  "Neuron 7",  "Neuron 8", ...
      "Neuron 9",  "Neuron 10");
xlabel("Time")
ylabel("Neuron activity")

figure;
surf(X,Y,X3)
title("Normalized Activity of Neurons in a Recurrent Neural Network with Signal" + ...
    " Function 3 Over Time - Surface Plot")
xlabel("Neuron Index")
ylabel("Time")
zlabel("Neuron activity")

%final normalized f3 vs index
figure;
plot(X3(101, :), "LineWidth",2) 
title("The Final Normalized Activity for Neurons in a Recurrent Neural Network with Signal Function 3")
xlabel("Neuron Index")
ylabel("Neuron Activity")
xlim([1 10])
% ylim([.001 .002]) 

%normalized f4 vs time - each with 10 curves
figure;
plot(times,X4,'LineWidth',2)
title("Normalized Activity of Neurons in a Recurrent Neural Network with Signal" + ...
    " Function 4 Over Time")
legend("Neuron 1", "Neuron 2", "Neuron 3",  "Neuron 4", ...
     "Neuron 5",  "Neuron 6",  "Neuron 7",  "Neuron 8", ...
      "Neuron 9",  "Neuron 10");
xlabel("Time")
ylabel("Neuron activity")

figure;
surf(X,Y,X4)
title("Normalized Activity of Neurons in a Recurrent Neural Network with Signal" + ...
    " Function 4 Over Time - Surface Plot")
xlabel("Neuron Index")
ylabel("Time")
zlabel("Neuron activity")

%final normalized f4 vs index
figure;
plot(X4(101, :), "LineWidth",2) 
title("The Final Normalized Activity for Neurons in a Recurrent Neural Network with Signal Function 4")
xlabel("Neuron Index")
ylabel("Neuron Activity")
xlim([1 10])




