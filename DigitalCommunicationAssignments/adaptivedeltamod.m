clc;

%% Adaptive delta modulation on sinusoidal signal

% Sampling parameters
Fs = 8000;                  
dt = 1/Fs;                   
StopTime = 0.25;             
t = (0:dt:StopTime-dt)';  

% Sine wave:
Fc = 60;                     
x = sin(2*pi*Fc*t);
x = 10*x(1:133);

subplot(4,1,1);
plot(x);
xlabel('samples');
title('Adaptive Delta Modulated Sinusoidal Signal');
hold on;

% initial output of staircase (approximation)
xr = zeros(1,length(x)); 
st = zeros(1,length(x));
st(1) = 0.8;
err = zeros(1,length(x));
err(1)=-1;

for i=1:length(x)-1
    
    if xr(i)<=x(i)
        err(i+1)=1;
        st(i+1) = abs(st(i))*err(i+1) + st(1)*err(i);
    else
        err(i+1)=-1;
        st(i+1) = abs(st(i))*err(i+1) + st(1)*err(i);
    end
    xr(i+1)=xr(i) + st(i+1);
end
stairs(xr);
hold off

%demodulation

for i=1:st
    if st(i)>xr(i)
        err(i+1)=1;
        st(i+1) = abs(st(i))*err(i+1) + st(1)*err(i);
        xr(i+1)=xr(i) + st(i+1);
    else
        err(i+1)=-1;
        st(i+1) = abs(st(i))*err(i+1) + st(1)*err(i);
        xr(i+1)=xr(i) + st(i+1);
    end
end

hd = lowpass(xr,2*pi,100);
subplot(4,1,2)
plot(hd,'c');
xlabel('samples');
title('Demodulation of ADM Sinusoidal Signal');

%% Adaptive delta modulation on audio file

[signal,fs]=audioread('band.wav');
x = 100*signal(400:520);

subplot(4,1,3);
plot(x);
xlabel('samples');
title('Adaptive Delta Modulated Audio Signal');
hold on;

% initial output of staircase (approximation)
xr = zeros(1,length(x)); 
st = zeros(1,length(x));
st(1) = 0.8;
err = zeros(1,length(x));
err(1)=-1;

for i=1:length(x)-1
    
    if xr(i)<=x(i)
        err(i+1)=1;
        st(i+1) = abs(st(i))*err(i+1) + st(1)*err(i);
    else
        err(i+1)=-1;
        st(i+1) = abs(st(i))*err(i+1) + st(1)*err(i);
    end
    xr(i+1)=xr(i) + st(i+1);
end
stairs(xr);
hold off

%demodulation

for i=1:st
    if st(i)>xr(i)
        err(i+1)=1;
        st(i+1) = abs(st(i))*err(i+1) + st(1)*err(i);
        xr(i+1)=xr(i) + st(i+1);
    else
        err(i+1)=-1;
        st(i+1) = abs(st(i))*err(i+1) + st(1)*err(i);
        xr(i+1)=xr(i) + st(i+1);
    end
end

hd = lowpass(xr,2*pi,100);
subplot(4,1,4)
plot(hd,'c');
xlabel('samples');
title('Demodulation of ADM Audio Signal');