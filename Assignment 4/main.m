clc;
clear all;
close all;

%% Utilities

snr_db = 0:10;
snr_lin=10.^(snr_db./10);
N = 10.^6;
b = randi([0 1], 1, N); %Txd Bits

%% Theoretical
bpsk_theo = zeros(1,length(snr_db));
qpsk_theo = zeros(1,length(snr_db));

for i=1:length(snr_lin)
    bpsk_theo(i) = qfunc(sqrt(2*snr_lin(i)));
    qpsk_theo(i) = qfunc(sqrt(2*snr_lin(i)));
end

%% Function Calls and Plotting

[bpsk_ber_prac, bpsk_ser_prac] = bpsk(b, snr_db);
[qpsk_ber_prac, qpsk_ser_prac] = qpsk(b, snr_db);

%% All Plotting

semilogy(snr_db,bpsk_theo,'-x');
hold on;
semilogy(snr_db,bpsk_ber_prac,'-o');
hold on;
semilogy(snr_db,bpsk_ser_prac,'-square');
hold on;
semilogy(snr_db,qpsk_theo,'-diamond');
hold on;
semilogy(snr_db,qpsk_ber_prac,'-*');
hold on;
semilogy(snr_db,qpsk_ser_prac,'-');
grid on;
title('BER vs log(SNR)');
ylabel('log(BER)')
xlabel('SNR')
legend('Theoretical BER BPSK', 'Simulated BER BPSK', 'Simulated SER BPSK', 'Theoretical BER QPSK', 'Simulated BER QPSK', 'Simulated SER QPSK')

%% BPSK Plotting

semilogy(snr_db,bpsk_theo,'-x');
hold on;
semilogy(snr_db,bpsk_ber_prac,'-o');
hold on;
semilogy(snr_db,bpsk_ser_prac,'-square');
grid on;
title('BER vs log(SNR) - BPSK');
ylabel('log(BER)')
xlabel('SNR')
legend('Theoretical BER', 'BER BPSK', 'SER BPSK')

%% QPSK Plotting

semilogy(snr_db,qpsk_theo,'-diamond');
hold on;
semilogy(snr_db,qpsk_ber_prac,'-*');
hold on;
semilogy(snr_db,qpsk_ser_prac,'-');
grid on;
title('BER vs log(SNR) - QPSK');
ylabel('log(BER)')
xlabel('SNR')
legend('Theoretical BER', 'BER QPSK', 'SER QPSK')

%% Functions

function [bpsk_ber, bpsk_ser] = bpsk(b, snr_db)
    snr_lin=10.^(snr_db./10);
    N = length(b);
    s = zeros(1, N);
    b_hat = zeros(1, N);
    s_hat = zeros(1, N);
    y = zeros(1, N);
    bpsk_ber = zeros(1,length(snr_db));
    bpsk_ser = zeros(1,length(snr_db));
    for j=1:length(snr_db)
       Es = 1;
       sd = sqrt(Es/(2*snr_lin(j)));
       noise = sd*(randn(1,N) + 1i*randn(1, N));
       for jj=1:length(b)
           if b(jj)==0
               s(jj)=1;
           else
               s(jj)=-1;
           end
       end
       y=s/sqrt(Es)+noise;
       for jj=1:length(b)
           if abs(y(jj)-1) > abs(y(jj)+1)
               s_hat(jj)=-1;
               b_hat(jj)=1;
           else
               s_hat(jj)=1;
               b_hat(jj)=0;
           end
       end
       biterror=xor(b_hat,b);
       symerror=abs(s - s_hat)/2;
       bpsk_ber(j)=sum(biterror)/N;
       bpsk_ser(j)=sum(symerror)/N;
    end
end


function [qpsk_ber, qpsk_ser] = qpsk(b, snr_db)
    snr_lin=10.^(snr_db./10);
    N = length(b);
    decs = zeros(1, N/2);
    s = zeros(1, N/2);
    y = zeros(1, N/2);
    b_hat = zeros(1, N);
    s_hat = zeros(1, N/2);
    mappings = zeros(1, N/2);
    qpsk_ber = zeros(1,length(snr_db));
    qpsk_ser = zeros(1,length(snr_db));
    odd=b(1:2:end);
    even=b(2:2:end);
    phi=[1+1i, -1+1i, 1-1i, -1-1i]/sqrt(2);
    maps = [0, 1, 2, 3];
    Es=1;
    
    for j=1:length(snr_db)
       sigma_2 = Es/(2*snr_lin(j));
       noise = sqrt(sigma_2/2)*(randn(1,N/2) + 1i*randn(1, N/2));
       % Bits to symbol mapping
       for jj=1:length(odd)
           decs(jj)=bin2dec(strcat(num2str(odd(jj)), num2str(even(jj))));
           s(jj)=phi(decs(jj)+1);
       end
       % Passing symbols through AWGN Channel
       y=s+noise;
       for jj=1:length(odd)
           d = abs(y(jj)-phi).^2;
           I=find(d==min(d));
           s_hat(jj)=phi(I);
           mappings(jj)=maps(I);
       end
       binaries = de2bi(mappings, 2, 'left-msb');
       b_hat= reshape(binaries.',1,[]);
       symerror = nnz(abs(s - s_hat));
       qpsk_ser(j) = 2*sum(symerror)/N;
       errors = xor(b_hat, b);
       qpsk_ber(j) = sum(errors)/(N);
    end
end