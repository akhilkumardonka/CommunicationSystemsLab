clc;
clear all;

%% Generating Practical BERs

SNR_db = 0:10;
SNR_lin=10.^(SNR_db./10);
energy = 1;
noises = energy./SNR_lin;
NA = 1000;
pErrTheo = zeros(1,length(SNR_db));
pErrSim50 = zeros(1,length(SNR_db));

% Theoretical
for i=1:length(SNR_lin)
    pErrTheo(i) = qfunc(sqrt(2*energy/noises(i)));
end

% Generated Enough Number of Blocks
NF = 1000;
bits = randi([0 1], NA, NF);
symbols = zeros(NA, NF);
r = zeros(NA, NF);
a_hat = zeros(NA, NF);

% Simulated
for k=1:length(noises)
    N0 = noises(k);
    p = zeros(1, NF);
    for i = 1:NF
        symbols(:, i) = (1 - 2*bits(:, i))*sqrt(energy);
        noise  = sqrt(N0/2)*randn(NA,1);
        r(:,i) = symbols(:, i) + noise;
    
        % Decision
        for j=1:NA
            if r(j, i)>=0
                a_hat(j, i)=0;
            else
                a_hat(j, i)=1;
            end
        end
        
        % Blockwise errors
        error = xor(bits(:, i), a_hat(:, i));
        p(i) = sum(error);
    end
    pErrSim50(k) = sum(p)/(NA*NF);
end

semilogy(SNR_db,pErrTheo);
hold on;
semilogy(SNR_db,pErrSim50);
title('BER vs log(SNR) - Theoretical and Practical');
ylabel('log(BER)')
xlabel('SNR')
legend('Theoretical BPSK', 'BPSK Simulated')