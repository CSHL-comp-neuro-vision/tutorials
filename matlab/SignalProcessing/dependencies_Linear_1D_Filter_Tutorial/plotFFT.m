function plotfft(t,y)


Y = complex2real(fft(y),t);

clf
subplot(1,2,1)
plot(t,y,'b-');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(1,2,2)
stem(Y.freq,Y.amp,'fill');
set(gca,'XLim',[0,max(Y.freq)]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');




