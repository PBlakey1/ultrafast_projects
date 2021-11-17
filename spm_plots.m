%Fourier Transform Parameters
N = 5001;
time = linspace(-20*25e-9,20*25e-9,N);
Ts = mean(diff(time));
Fs = 1/Ts;
N2 = 2^nextpow2(N);
Ts = 1/Fs;
F = Fs/N2;
Nyq = Fs/2;
freq = -Nyq:F:Nyq-F;
freq2 = linspace(-1,1,fix(numel(time)))*Nyq;




%Physical Parameters
c = 2.998e8;
lambda0 = 1550e-9;
w0 = 2*pi*c/lambda0;
n0 = 1.5;
n2 = 3e-20;
%epsilon0 = 8.854e-12;
gamma = w0*n2/c;
%gamma = 0.9e-10;
T0 = 6e-9; %ps
a0 =sqrt(1e13);
L = 10;

%Pulse Propagation
a_t = a0*sech(time/T0);
delta_phi = gamma.*a_t.*conj(a_t)*L;
num_peaks = max(delta_phi)/pi + 1;
aL_t = a_t.*exp(-1i.*delta_phi);
I0_t = a_t.*conj(a_t);
IL_t = aL_t.*conj(aL_t);

%Instantanious Frequencies
inst_freq_0 = w0*ones(1,length(a_t));
inst_freq_L = -gradient(delta_phi(:))./gradient(time(:));

%Fourier Transforms
a_w = fftshift(fft(a_t)/numel(a_t));
aL_w = fftshift(fft(aL_t)/numel(aL_t));
I0_w = abs(a_w)*2;
IL_w = aL_w.*conj(aL_w);

%Plots
subplot(4,2,1)
plot(time, I0_t,'-r');
xlabel('Time(s)');
ylabel('Intensity (W/m^2)')
subplot(4,2,2)
plot(time, IL_t,'-r');
xlabel('Time(s)');
ylabel('Intensity (W/m^2)')
subplot(4,2,3)
plot(time, inst_freq_0, '-');
subplot(4,2,4)
plot(time, inst_freq_L, '-');

subplot(4,2,5)
plot(freq2, I0_w, '-');
xlim([-1e9,1e9])
subplot(4,2,6)
plot(freq2, IL_w, '-');
xlim([-1e9,1e9])

subplot(4,2,7)
plot(time, delta_phi, '-');
