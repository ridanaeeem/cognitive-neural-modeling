%%% part a 
close all
A = 0.1;
B = 1;
x = zeros(1,10);
I = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];
time = [0,1,2,3,4,5,6,7,8,9,10];
index2 = [0,1,2,3,4,5,6,7,8,9];

%%additive

dxdt = [];
%get x values for each time through Euler's method
for i=1:9
    sumExceptI = sum(I) - I(i);
    disp(B-sumExceptI)
    dxdt(i) = (-A * x(i)) + (B * I(i)) - sumExceptI;
    x(i+1) = x(i) + (dxdt(i));
end
disp("additive x")
disp(x)
%plot x
%FIGURE 1
figure;
plot(index2,x,'LineWidth',2)
hold on
plot(index2,I,'LineWidth',2)
title("Effect of Stimulus on Activity of a Neuron in a Feedforward Additive Network Over Time")
xlabel("Time")
ylabel("Activity")
legend("Activity","Input")
%normalize + plot X
X = normalized(x);
disp("normalized additive")
disp(X)
%FIGURE 2
figure;
plot(index2,X,'LineWidth',2)
hold on
plot(index2,I,'LineWidth',2)
title("Normalized Activity of Additive Feedforward Network in Response to Stimulus Over Time")
xlabel("Time")
ylabel("Activity")
legend("Activity","Input")

%%shunting
x = zeros(1,10);
dxdt=[];
for i=1:9
    sumExceptI = sum(I) - I(i);
    dxdt(i) = (-A * x(i)) + ((B - x(i)) * I(i)) - (x(i) * sumExceptI);
    x(i+1) = x(i) + (dxdt(i));
end
disp("shunting x")
disp(x)
%plot x
%FIGURE 3
figure;
plot(index2,x,'LineWidth',2)
hold on
plot(index2,I,'LineWidth',2)
title("Effect of Stimulus on Activity of a Neuron in a Feedforward Shunting Network Over Time")
xlabel("Time")
ylabel("Activity")
legend("Activity","Input")
%normalize + plot X
X = normalized(x);
disp("normalized shunting")
disp(X)
%FIGURE 4
figure;
plot(index2,X,'LineWidth',2)
hold on
plot(index2,I,'LineWidth',2)
xlabel("Time")
ylabel("Activity")
title("Normalized Activity of Shunting Feedforward Network in Response to Stimulus Over Time")
legend("Activity","Input")
%%
%%%part c
%close all
A = 0.1;
B = 1;
x = zeros(1,10);
dxdt = [];
I = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1];
time = [0,1,2,3,4,5,6,7,8,9,10];
index2 = [0,1,2,3,4,5,6,7,8,9];

%get x values for each time through Euler's method
for i=1:9
    sumExceptI = sum(I) - I(i);
    disp(B-sumExceptI)
    dxdt(i) = (-A * x(i)) + ((B - x(i)) * I(i)) - sumExceptI;
    x(i+1) = x(i) + (dxdt(i));
end

%plot x
%FIGURE 5
figure;
plot(index2,x,'LineWidth',2)
hold on
plot(index2,I,'LineWidth',2)
title("Effect of Stimulus on Activity of a Neuron in a Modified Feedforward Additive Network")
xlabel("Time")
ylabel("Activity")
legend("Activity","Input")

%normalize + plot X
X = normalized(x);
disp("normalized shunting")
disp(X)
%FIGURE 6
figure;
plot(index2,X,'LineWidth',2)
hold on
plot(index2,I,'LineWidth',2)
xlabel("Time")
ylabel("Activity")
title("Normalized Activity of Shunting Feedforward Network in Response to Stimulus Over Time")
legend("Activity","Input")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%part d
close all
A = 0.1;
B = 1;
x = zeros(1,10);
dxdt = [];
time = [0,1,2,3,4,5,6,7,8,9,10];
index2 = [0,1,2,3,4,5,6,7,8,9];

I1 = [0.1,0.1,0.1,0.1,0.1,0.8,0.8,0.8,0.8,0.8];
I2 = [0.1,0.1,0.1,0.1,0.8,0.8,0.1,0.1,0.1,0.1];
I3 = [0.1,0.1,0.8,0.8,0.8,0.8,0.8,0.8,0.1,0.1];
I4 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];

for i=1:9
    sumExceptI = sum(I4) - I4(i);
    %disp(B-sumExceptI)
    dxdt(i) = (-A * x(i)) + ((B - x(i)) * I4(i)) - (x(i) * sumExceptI);
    x(i+1) = x(i) + (dxdt(i));
end

disp("compare d to this old one without coeffs")
disp(x)

%plot x
% figure;
% plot(index2,x,'LineWidth',2)
% hold on
% plot(index2,I1,'LineWidth',2)
% title("Shunting Equation 1 with Input A")
% xlabel("Time")
% ylabel("Activity")
% legend("Activity","Input")

%normalize + plot X
X = normalized(x);
disp("normalized shunting")
disp(X)
%FIGURE 
index3 = [5,6,7,8,9,10,11,12,13,14];
figure;
plot(index3,X,'LineWidth',2)
hold on
xlabel("Time")
ylabel("Activity")

%modified
I1 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8];
I2 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.8,0.8,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1];
I3 = [0.1,0.1,0.1,0.1,0.1,0.1,0.8,0.8,0.8,0.8,0.8,0.8,0.1,0.1,0.1,0.1,0.1,0.1];
I4 = [0.1,0.1,0.1,0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.0,1.0,1.0,1.0];

x = zeros(1,10);
dxdt = [];
A = 0.1;
B = 10;

%new shunting - I1
for i=1+4:9+4
    k = i - 4;
    t = i + 4;
    prodC = 0;
    prodE = 0;
    sumC = 0;
    sumE = 0;
    for j=k:t
        text = [' this is number ',num2str(j), ' out of ', num2str(t -k)];
        disp(text)
        curC = exp((-(j-i)^2)/4)
        curE = 0.5 * exp((-(j-i)^2)/16)
        prodC = curC * I4(j)
        prodE = curE * I4(j)
        sumC = sumC + prodC
        sumE = sumE + prodE
    end
    dxdt(i) = (-A * x(i)) + ((B - x(i)) * prodC) - (x(i) * prodE)
    x(i+1) = x(i) + (dxdt(i))
end
%plot x
% figure;
% plot(x,'LineWidth',2)
% hold on
% plot(I1,'LineWidth',2)
% title("shunt2 I1")
% xlabel("Time")
% ylabel("Activity")
% legend("Activity","Input")

%normalize + plot X
X = normalized(x);
disp("normalized shunting")
disp(X)
%FIGURE 6
% figure;
plot(X,'LineWidth',2)
hold on
plot(I4,'LineWidth',2)
xlabel("Time")
ylabel("Activity")
title("Normalized Activity of Two Feedforward Shunting Networks in Response to Stimulus D Over Time")
legend("Equation 3","Equation 4","Input D")


%%%%%%%%%%%%%%%
%%

%new shunting - I2
for i=1+4:9+4
    k = i - 4;
    t = i + 4;
    prodC = 0;
    prodE = 0;
    sumC = 0;
    sumE = 0;
    for j=k:t
        text = [' this is number ',num2str(j), ' out of ', num2str(t -k)];
        disp(text)
        curC = exp((-(j-i)^2)/4)
        curE = 0.5 * exp((-(j-i)^2)/16)
        prodC = curC * I2(j)
        prodE = curE * I2(j)
        sumC = sumC + prodC
        sumE = sumE + prodE
    end
    dxdt(i) = (-A * x(i)) + ((B - x(i)) * prodC) - (x(i) * prodE)
    x(i+1) = x(i) + (dxdt(i))
end
%plot x
figure;
plot(x,'LineWidth',2)
hold on
plot(I2,'LineWidth',2)
title("shunt I2")
xlabel("Time")
ylabel("Activity")
legend("Activity","Input")


x = zeros(1,10);
dxdt = [];

%new shunting - I1
for i=1+4:9+4
    k = i - 4;
    t = i + 4;
    prodC = 0;
    prodE = 0;
    sumC = 0;
    sumE = 0;
    for j=k:t
        text = [' this is number ',num2str(j), ' out of ', num2str(t -k)];
        disp(text)
        curC = exp((-(j-i)^2)/4)
        curE = 0.5 * exp((-(j-i)^2)/16)
        prodC = curC * I3(j)
        prodE = curE * I3(j)
        sumC = sumC + prodC
        sumE = sumE + prodE
    end
    dxdt(i) = (-A * x(i)) + ((B - x(i)) * prodC) - (x(i) * prodE)
    x(i+1) = x(i) + (dxdt(i))
end
%plot x
figure;
plot(x,'LineWidth',2)
hold on
plot(I3,'LineWidth',2)
title("shunt I3")
xlabel("Time")
ylabel("Activity")
legend("Activity","Input")


x = zeros(1,10);
dxdt = [];

%new shunting - I1
for i=1+4:9+4
    k = i - 4;
    t = i + 4;
    prodC = 0;
    prodE = 0;
    sumC = 0;
    sumE = 0;
    for j=k:t
        text = [' this is number ',num2str(j), ' out of ', num2str(t -k)];
        disp(text)
        curC = exp((-(j-i)^2)/4)
        curE = 0.5 * exp((-(j-i)^2)/16)
        prodC = curC * I4(j)
        prodE = curE * I4(j)
        sumC = sumC + prodC
        sumE = sumE + prodE
    end
    dxdt(i) = (-A * x(i)) + ((B - x(i)) * prodC) - (x(i) * prodE)
    x(i+1) = x(i) + (dxdt(i))
end
%plot x
figure;
plot(x,'LineWidth',2)
hold on
plot(I4,'LineWidth',2)
title("shunt I4")
xlabel("Time")
ylabel("Activity")
legend("Activity","Input")

%%
%%C AND E
y=[];
for i=0.1:0.1:10
    pos=int16(i/.1);
    y(pos) = exp(-(i^2) / 4);
end
figure;
plot(y,'LineWidth',2)

z=[];
for i=0.1:0.1:10
    pos=int16(i/.1);
    z(pos) = 0.5 * exp(-(i^2) / 16);
end
hold on
plot(z,'LineWidth',2)
title("Varying Distance-Dependent Receptive Field Coefficient Values for Changing |k-i|")
ylabel("Value")
xlabel("|k - i|")
legend("C","E")
%% defining functions
%additive for part a
function dxdt = additive(A, x, B, I)
    for i=1:10
        sumExceptI = sum(I) - I(i);
        dxdt(i) = (-A * x(i)) + (B * I(i)) - sumExceptI;
    end
end

%shunting for part a
function dxdt = shunting(A, x, B, I)
    for i=1:10
        sumExceptI = sum(I) - I(i);
        dxdt(i) = (-A * x(i)) + (I(i) * (B - x(i))) - (x(i) * sumExceptI);
    end
end

%normalized for part a
function X = normalized(x)
    for i=1:size(x,2)
        X(i) = x(i)/sum(x);
    end
end

%additive for part c
function dxdt = additive2(A, x, B, I)
    for i=1:10
        sumExceptI = sum(I) - I(i);
        dxdt(i) = (-A * x(i)) + ((B - x(i)) * I(i)) - sumExceptI;
    end
end

%shunting for part d
function dxdt = shunting2(A, x, B, I, C, E)
    for i=1:10
        sumExceptI = sum(I) - I(i);
        dxdt(i) = (-A * x(i)) + (I(i) * (B - x(i))) - (x(i) * sumExceptI);
    end
end