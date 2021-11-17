%Fourier Transform Parameters
N = 5001;
time = linspace(-25,25,N);
Ts = mean(diff(time));
Fs = 1/Ts;
N2 = 2^nextpow2(N);
Ts = 1/Fs;
F = Fs/N2;
Nyq = Fs/2;
freq = -Nyq:F:Nyq-F;




%Physical Parameters
c = 2.998e8;
lambda0 = 1550e-6;
w0 = 2*pi*c/lambda0;
n0 = 1.5;
n2 = 3e-20;
epsilon0 = 8.854e-12;
gamma = w0*epsilon0*n0*n2;
gamma = 0.9e-4;
T0 = 6; %ps
a0 =sqrt(1e13);
L = 10;

%Pulse Propagation
delta_phi = gamma.*a_t.*conj(a_t)*L;

a_t = a0*sech(time/T0);
aL_t = a_t.*exp(-1i.*delta_phi);
I0_t = a_t.*conj(a_t);
IL_t = aL_t.*conj(aL_t);

%Instantanious Frequencies
inst_freq_0 = w0*ones(1,length(a_t));
inst_freq_L = - gradient(delta_phi);

%Fourier Transforms
a_w = fftshift(fft(a_t,N2));
aL_w = fftshift(fft(aL_t,N2));
I0_w = a_w.*conj(a_w);
IL_w = aL_w.*conj(aL_w);

%Plots
subplot(3,2,1)
plot(time, I0_t,'-r');
xlabel('Time(ps)');
ylabel('Intensity (W/m^2)')
subplot(3,2,2)
plot(time, IL_t,'-r');
xlabel('Time(ps)');
ylabel('Intensity (W/m^2)')
subplot(3,2,3)
plot(time, inst_freq_0, '-');
subplot(3,2,4)
plot(time, inst_freq_L, '-');

subplot(3,2,5)
plot(freq, I0_w, '-');
subplot(3,2,6)
plot(freq, IL_w, '-');