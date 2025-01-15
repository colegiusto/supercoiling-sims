function dxdt = ODE(~,x)
    % states are in groups of 3, with 
    % 1,2,3 being related to topisomerase I,
    % 4,5,6 is gyrase
    % 7,8,9 is Fis
    % 10, 11, 12 is cspA

    %The first variable is transcriptional state, the second is protein
    %state, the third is supercoiling state at the promoter
    
    N = 4; % number of genes in the network

    alpha = ones(N,1); % gains for mRNA transcription each gene in order
    dm = ones(N,1); % rates of mRNA decay

    a_p = 1; % rate of protein transcription from mRNA
    dp = 1; % rate of protein decay

    a_sigma = 1; %rate of positive supercoiling from upstream genes
    
    
     
    
    s = ones(N,1); % Scale of rbfs

    k_sigma = @(sig, sig_opt)exp((sig-sig_opt)^2/s(1)); %supercoiling effect function

    k = 1; %Gene network edge strength

    sigma_star = [0 0 -1 -1]; % optimal supercoiling densities
    
    
    
    
    
    dxdt = zeros(length(x),1);
    
    dxdt(1) = alpha(1)*k_sigma(x(3), sigma_star(1))*1/(1+k*x(8))-dm(1)*x(1);                   %mRNA for topisomerase
    dxdt(4) = alpha(2)*k_sigma(x(6), sigma_star(2))*(k*x(11)/(1+k*x(8)+k*x(11)))-dm(2)*x(4);   %mRNA for Gyrase
    dxdt(7) = alpha(3)*k_sigma(x(9), sigma_star(3))*(k*x(11)/(1+k*x(11)))-dm(3)*x(7);          %mRNA for Fis
    dxdt(10) = alpha(4)*k_sigma(x(12), sigma_star(4))*(k*x(8)/(1+k*x(8)))-dm(4)*x(10);         %mRNA for cspA
    
    for i = 1:N
        j = 3*(N-1);
        dxdt(j+2) = a_p*x(j+1)-dp*x(j+2);
        dxdt(j+3) = -(x(2)+x(5))*x(j+3)-dxdt(j+1)-dm(i)*x(j+1)+a_sigma;
    end

end