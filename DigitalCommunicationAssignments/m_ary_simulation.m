clc;
clear all;
close all;

%% Generating Bits

M_arrs = [2, 4, 16, 64, 256];
SNR_db = 0:2:16;
SNR_lin = 10.^(SNR_db/10);
NA = 10080;
NF = 100;
data = randi([0 1], NF, NA);
ASK_theo = zeros(length(M_arrs), length(SNR_lin));
PSK_theo = zeros(length(M_arrs), length(SNR_lin));
QAM_theo = zeros(length(M_arrs), length(SNR_lin));

%% Theoretical Plots

for k=1:length(M_arrs)
    M = M_arrs(k);
    for i = 1:length(SNR_db)
        Eb_N0 = SNR_lin(i);
        ASK_theo(k, i) = ((M-1)/(M*log2(M)))*qfunc(sqrt((3*log2(M)/(M^2-1))*(Eb_N0)));
        PSK_theo(k, i) = (2/log2(M))*qfunc(sqrt(2*Eb_N0*log2(M))*sin(pi/M));
        QAM_theo(k, i) = 2*((sqrt(M)-1)/(sqrt(M)*log2(M)))*qfunc(sqrt((3*log2(M)*Eb_N0)/(2*(M-1))));
    end
end

%% Simulated ASK

ASK_sim = zeros(length(M_arrs), length(SNR_lin));

for m=1:length(M_arrs)
    M = M_arrs(m);
    k = log2(M); 
    alphabet = (-(M-1) : 2 : (M-1));
    Eavg = (1/M)*(sum(alphabet.^2)); 
    Eb_N0 = SNR_lin;
    Es_N0 = Eb_N0*k;
    
    St_blocks = zeros(NF, NA/k);

    % bit to symbol mapping
    for b=1:NF
        block = data(b,:);
        syms = bi2de(reshape(block,k,length(block)/k).','left-msb');
        symbols = [];
        for s=1:length(syms)
            symbols = [symbols alphabet(syms(s)+1)];
        end
        St_blocks(b,:)=symbols;
    end
    
    % Txn through AWGN
    
    St = St_blocks(1,:);
    St_norm = St/sqrt(Eavg);
    sigma2 = 1./(Es_N0);                
    N0 = sigma2;
    Sr_norm = zeros(length(SNR_db),NA/k);
    for p = 1:length(SNR_db)
        n = sqrt(sigma2(p)/2)*(randn(1,length(St))+1i*randn(1,length(St)));
        Sr_norm(p,:) = St_norm + n;
    end
    Sr = Sr_norm*sqrt(Eavg);               
    DM = [-(M-2):2:(M-2)];
    So_I = zeros(1,NA/k);
    So = zeros(length(SNR_db),NA/k);
    for i = 1:length(SNR_db)
        So_I(find(real(Sr(i,:)) < DM(1))) = alphabet(1);
        if (length(DM) > 1)
            for k = 2:length(DM)
                So_I(find((real(Sr(i,:)) >  DM(k-1))&(real(Sr(i,:)) < DM(k)))) = alphabet(k);
            end
        end
        So_I(find(real(Sr(i,:)) > DM(length(DM)))) = alphabet(length(alphabet));
        So(i,:) = So_I;
        
        ASK_sim(m,i) = symerr(St,So(i,:))/NA;
    end

end

%% MPSK Simulated Plots

PSK_sim = zeros(length(M_arrs), length(SNR_lin));

for m=1:length(M_arrs)
    M = M_arrs(m);
    k = log2(M);
    N = NA/k;
    alphabet = (0:M-1)*2*pi/M;
    sdB  = SNR_db + 10*log10(k);
    b = (0:M-1);
    map = bitxor(b,floor(b/2));
    [tt, ind] = sort(map);                                
    c = zeros(1,N);
    errors = zeros(length(SNR_db),1);
    for i = 1:length(SNR_db)
        for j=1:NF
            bits = data(j,:);
            bin2DecMatrix = ones(N,1)*(2.^((k-1):-1:0)) ;   
            shape = reshape(bits,k,N).';
            G= (sum(shape.*bin2DecMatrix,2)).';
            dec = ind(G+1)-1; 
            ph= dec*2*pi/M;
            d= exp(1i*ph);   
            s = d;
            n = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N)); 
            r = s + 10^(-sdB(i)/20)*n; 
            e = angle(r);     
            e(e<0) = e(e<0) + 2*pi;
            c = 2*pi/M*round(e/(2*pi/M));    
            c(c==2*pi) = 0;
            cd = round(c*M/(2*pi));        
            f = map(cd+1); 
            cb = dec2bin(f,k) ;      
            cb = cb.';
            cb = cb(1:end).';
            cb = str2num(cb).' ;
            Err(i,j) = size(find(bits- cb),2);
        end
        errors(i) = sum(Err(i,:))/NF;
    end
    PSK_sim(m,:)=errors/(N*k);
end

%% MQAM Simulated

QAM_sim = zeros(length(M_arrs), length(SNR_lin));

for marr=1:length(M_arrs)
    M = 16;
    k = log2(M);
    N = NA/k;
    Re = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
    Im = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
    k_QAM = 1/sqrt(10);
    sdB  = SNR_db + 10*log10(k);
    a = (0:k-1);
    map = bitxor(a,floor(a/2));
    [tt, ind] = sort(map);
    errors = zeros(length(SNR_db),1);
    c = zeros(1,N);
    for i = 1:length(SNR_db)
        for j=1:NF
            c = data(1,:);
            d = reshape(c,k,N).';
            bd = ones(N,1)*(2.^((k/2-1):-1:0));
            cRe =  d(:,(1:k/2));
            e = sum(cRe.*bd,2);
            f = bitxor(e,floor(e/2));
            cIm =  d(:,(k/2+1:k));
            g = sum(cIm.*bd,2);
            h = bitxor(g,floor(g/2));
            modRe = Re(f+1);
            modIm = Im(h+1);
            mod = modRe + 1i*modIm;
            s = k_QAM*mod;
            n = 1/sqrt(2)*[randn(1,N) + 1i*randn(1,N)];  
            r = s + 10^(-sdB(i)/20)*n; 
            r_re = real(r)/k_QAM; 
            r_im = imag(r)/k_QAM;
            m = 2*floor(r_re/2)+1;
            m(m>max(Re)) = max(Re);
            m(m<min(Re)) = min(Re);
            n= 2*floor(r_im/2)+1;
            n(n>max(Im)) = max(Im);
            n(n<min(Im)) = min(Im);
            oRe = ind(floor((m+4)/2+1))-1; 
            oIm = ind(floor((n+4)/2+1))-1;
            pRe = dec2bin(oRe,k/2);
            pIm = dec2bin(oIm,k/2);
            pRe = pRe.';
            pRe = pRe(1:end).';
            pRe = reshape(str2num(pRe).',k/2,N).';
            pIm = pIm.';
            pIm = pIm(1:end).';
            pIm = reshape(str2num(pIm).',k/2,N).';
            Err(i,j) = size(find([cRe- pRe]),1) + size(find([cIm - pIm]),1) ;
        end
        errors(i) = sum(Err(i,:))/NF;
    end
   QAM_sim(marr,:) = errors/(N*k);
end

%% Final Plots M = 2

semilogy(SNR_db,ASK_theo(1,:));
hold on;
semilogy(SNR_db,ASK_sim(1,:));
hold on;
semilogy(SNR_db,PSK_theo(1,:));
hold on;
semilogy(SNR_db,PSK_sim(1,:));
hold on;
semilogy(SNR_db,QAM_theo(1,:));
hold on;
semilogy(SNR_db,QAM_sim(1,:));
title('BER vs log(SNR) for M = 2');
ylabel('log(BER)')
xlabel('SNR')
legend('ASK Theoretical', 'ASK Simulated', 'PSK Theoretical', 'PSK Simulated','QAM Theoretical', 'QAM Simulated')

%% Final Plots M = 4

semilogy(SNR_db,ASK_theo(2,:));
hold on;
semilogy(SNR_db,ASK_sim(2,:));
hold on;
semilogy(SNR_db,PSK_theo(2,:));
hold on;
semilogy(SNR_db,PSK_sim(2,:));
hold on;
semilogy(SNR_db,QAM_theo(2,:));
hold on;
semilogy(SNR_db,QAM_sim(2,:));
title('BER vs log(SNR) for M = 4');
ylabel('log(BER)')
xlabel('SNR')
legend('ASK Theoretical', 'ASK Simulated', 'PSK Theoretical', 'PSK Simulated','QAM Theoretical', 'QAM Simulated')

%% Final Plots M = 16

semilogy(SNR_db,ASK_theo(3,:));
hold on;
semilogy(SNR_db,ASK_sim(3,:));
hold on;
semilogy(SNR_db,PSK_theo(3,:));
hold on;
semilogy(SNR_db,PSK_sim(3,:));
hold on;
semilogy(SNR_db,QAM_theo(3,:));
hold on;
semilogy(SNR_db,QAM_sim(3,:));
title('BER vs log(SNR) for M = 16');
ylabel('log(BER)')
xlabel('SNR')
legend('ASK Theoretical', 'ASK Simulated', 'PSK Theoretical', 'PSK Simulated','QAM Theoretical', 'QAM Simulated')

%% Final Plots M = 64

semilogy(SNR_db,ASK_theo(4,:));
hold on;
semilogy(SNR_db,ASK_sim(4,:));
hold on;
semilogy(SNR_db,PSK_theo(4,:));
hold on;
semilogy(SNR_db,PSK_sim(4,:));
hold on;
semilogy(SNR_db,QAM_theo(4,:));
hold on;
semilogy(SNR_db,QAM_sim(4,:));
title('BER vs log(SNR) for M = 64');
ylabel('log(BER)')
xlabel('SNR')
legend('ASK Theoretical', 'ASK Simulated', 'PSK Theoretical', 'PSK Simulated','QAM Theoretical', 'QAM Simulated')

%% Final Plots M = 256

semilogy(SNR_db,ASK_theo(5,:));
hold on;
semilogy(SNR_db,ASK_sim(5,:));
hold on;
semilogy(SNR_db,PSK_theo(5,:));
hold on;
semilogy(SNR_db,PSK_sim(5,:));
hold on;
semilogy(SNR_db,QAM_theo(5,:));
hold on;
semilogy(SNR_db,QAM_sim(5,:));
title('BER vs log(SNR) for M = 256');
ylabel('log(BER)')
xlabel('SNR')
legend('ASK Theoretical', 'ASK Simulated', 'PSK Theoretical', 'PSK Simulated','QAM Theoretical', 'QAM Simulated')