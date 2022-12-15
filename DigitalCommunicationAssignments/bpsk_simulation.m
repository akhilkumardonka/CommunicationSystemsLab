clc;
clear all;

%% Input Data Sequence

NA = 1000;
bits = randi([0 1], NA, 1);
symbols = zeros(NA, 1);

%% Symbol Mapper

energy = 1;

for n = 1:NA
    symbols(n) = (1 - 2*bits(n))*sqrt(energy);
end

%% AWGN Channel

N0 = 1;
noise  = sqrt(N0/2)*randn(NA,1);
r = symbols + noise;

%% Decision Device

a_hat = zeros(NA, 1);

for i=1:NA
    if r(i)>=0
        a_hat(i)=0;
    else
        a_hat(i)=1;
    end
end

%% Data sink

error = xor(bits,a_hat);
p_hat = sum(error)/NA;

%% Splitting to block level

NA = 1000;
NF = 100;

bits = randi([0 1], NA, NF);
symbols = zeros(NA, NF);
r = zeros(NA, NF);
p = zeros(1, NF);
a_hat = zeros(NA, NF);

energy = 1;
N0 = 1;

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
    p(1,i) = sum(error);
end

p_hat = sum(p)/(NA*NF);

%% Generating Practical BERs

SNR_db = 0:10;
SNR_lin=10.^(SNR_db./10);
energy = 1;
noises = energy./SNR_lin;
NA = 1000;
NF = 100;
bits = randi([0 1], NA, NF);
symbols = zeros(NA, NF);
r = zeros(NA, NF);
a_hat = zeros(NA, NF);
pErrSim = zeros(1,length(SNR_db));
pErrTheo = zeros(1,length(SNR_db));
pErrSim50 = zeros(1,length(SNR_db));

% Theoretical
for i=1:length(SNR_lin)
    pErrTheo(i) = qfunc(sqrt(2*energy/noises(i)));
end

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
    pErrSim(k) = sum(p)/(NA*NF);
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
semilogy(SNR_db,pErrSim);
hold on;
semilogy(SNR_db,pErrSim50);
title('BER vs log(SNR) - Theoretical and Practical');
ylabel('log(BER)')
xlabel('SNR')
legend('Theoretical', 'Practical 100 Blocks', 'Practical 1000 Blocks')