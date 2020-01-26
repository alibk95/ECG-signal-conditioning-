%%
%load the real ecg signal from database physionet
load('100m.mat');
signal = val(1,:);
plot(signal);
title('original signal');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulate an artificial ECG signal
h = 70;
p = 2.5;
t1=10;
x1 = p*ecg(6000).';
y1 = sgolayfilt(kron(ones(1,19),x1),0,25);
m = h / 6;
n1=75000;
n = 1:n1;
del = round(6000);
mhb = y1(n + del); 
t = linspace(0,t1,n1);
plot(t,mhb, 'g');
q=p+0.5;
axis([0 t1 -q q]); 
grid;
%%
%artificial noise
L = 75000; 
mu = 0;
sigma = 2; 
X = sigma + randn(L, 1) + mu;
figure(1)
plot(X);
new  = X + mhb;
%%
varn = 1:0.01:10;                            % Variance
y = 5*sin(4*pi*varn/5);                     % Signal
yn = y + sqrt(varn).*randn(size(varn))+1;   % Signal + Noise
figure(1)
plot(varn, yn)
grid
%%
%add noise 
figure(1)
y = awgn(mhb,20,'measured');
plot(y,'b')
grid;
xlabel('time [sec]');
ylabel('voltage [mV]');
title('ECG signal + Noise');
%%
%Add baseline drift
t = linspace(0,20*pi,75000);
sn = 0.5*cos(t./30);
figure(1)
subplot(2,1,1)
plot(t, mhb)
title('Raw ECG Signal')
grid
subplot(2,1,2)
noisy_signal = y + sn;
plot(t, noisy_signal)
title('Noisy ECG Signal')
grid
%%
B1 = [0 1 5 1 0];
B2 = [1 1 1 1 1];
var1 = imdilate(y, B1);
var2 = imerode(y, B1);

var11 = imerode(var1, B2);
var22  =imdilate(var2, B2);

ns_sig  = (var11 + var22)/2;
figure(2)
plot(ns_sig)
%%
% Baseline correction 
BO = ones(1, 1500);
BC = ones(1,2225);
opened_signal = imopen(noisy_signal, BO);
figure(1)
hold on
plot(opened_signal, 'b', 'LineWidth', 5)
plot(noisy_signal, 'g', 'LineWidth', 0.5)
legend('Opened Signal','Noisy ECG Signal');
ylim([-3 3]);
closed_signal = imclose(opened_signal, BC);
figure(2)
plot(closed_signal, 'b', 'LineWidth', 3)
legend('Baseline');
ylim([-3 5]);
baseline_corrected_signal = noisy_signal - closed_signal;
figure(3)
plot(baseline_corrected_signal);
title('Baseline corrected signal');
xlabel('time [sec]');
ylabel('voltage [mV]');
%again checking the baseline drift of the corrected signal 
% r = imopen(baseline_corrected_signal,BO);
% e = imclose(r, BC);
%figure
%plot(e, 'b', 'LineWidth', 3);
%%
%Noise suppression 
ns_sig = baseline_corrected_signal;
S1 = [0 1 5 1 0];
S2 = [0 0 0.0001 0 0];
iter = 100;
ns = zeros(iter,75000);
dff = zeros(iter-1, 75000);
ns(1, :) = ns_sig;
for i = 2 : iter
ar1 = imdilate(ns_sig, S1);
ar2 = imerode(ns_sig, S1);
ar11 = imerode(ar1, S2);
ar22  =imdilate(ar2, S2);
ns_sig  = (ar11 + ar22)/2;
ns(i,:) = ns_sig;
dff(i-1,:) = ns(i,:) - ns(i-1,:);
mean_error(i-1) = mse(ns_sig , mhb);
end
final = ns_sig;
figure(1)
plot(baseline_corrected_signal)
figure(2)
plot(final)
figure(3)
plot(mean_error)
xlabel('# Iteration');
ylabel('Mean square error (MSE)');
figure(4)
plot(mean(dff.^2, 2))
xlabel('# Iteration');
ylabel('Error');
min_mean_error = min(mean_error);
min_x_mean_error = find(mean_error == min_mean_error);
%%
%evaluation of the proposed method
%SDR smaller the better(signal distortion ratio)

%BCR Baseline correction ratio, the bigger the better 
BCR = sum(abs(closed_signal)) / sum(abs(baseline_corrected_signal));
%NSR Nosie suppression ratio, the bigger the better 
%%
%erosion
for i = abs((m - 1) / 2) : n - abs((m + 1) / 2) - 1
    for j = 1 : m 
        d(i) = min(d(i - 1) , ECGsignal(i - abs((m - 1) / 2) + j) - B(j));
    end
end
figure()
plot(d)
%%
%dilation
for i = abs((m - 1) / 2) : n - abs((m + 1) / 2) - 1
    for j = 1 : m 
        d(i) = max(d(i - 1) , ECGsignal(i - abs((m - 1) / 2) + j) - B(j));
    end
end
figure()
plot(d)
