function [CCent] = CausalCent(CSDs,PEB)

% Tahereh S. Zarghami (tszarghami{at}gmail)
% See LICENSE file

    n = sqrt(numel(PEB.Ep)); % #nodes
    Ns = numel(CSDs); % #subj
    Ff = PEB.F; % free energy of full model
    
    for node = 1:n 
    
        for subj = 1:Ns 
            
            % Full model info
            pC = CSDs{subj}.M.pC;
            pE = CSDs{subj}.M.pE;
            A = pE.A;
    
            % Specify indices of the connections to shrink
            shrink = zeros(size(A));
            shrink(:,node) = 1; % efferents = outgoing
            shrink(node,:) = 0; % afferents = incoming
            shrink_ind = find(shrink); % linear index
            noshrink_ind = find(shrink == 0);
        
            % Shrink selected connections to zero (with precise null priors)
            A(shrink_ind) = 0;
            rE = pE; rE.A = A;
            pC(shrink_ind,shrink_ind) = 0;
            pC(shrink_ind,:) = 0;
            pC(:,shrink_ind) = 0;
            rC = pC;
	        
	        % BMR on CSDs 
            CSD = CSDs{subj};
            CSDR = spm_dcm_reduce(CSD,rE,rC); % reduced CSD
            CSDRs {node,subj} = CSDR;
                
        end
    
        csdrs = CSDRs (node,:);
        
        % PEB on reduced CSDs
        PEB_R(node) = spm_dcm_peb(csdrs(:));
        Fr (node) = PEB_R(node).F; % free energy of reduced model
    
        % dF
        dF_PEB (node) =  Ff - Fr (node);
        
        % KLD
        p.mu = full(PEB.Ep);
        p.sigma = PEB.Cp;
        EpR = zeros(size(p.mu)); EpR(noshrink_ind) = full(PEB_R(node).Ep);
        q.mu = EpR;
        CpR = zeros(size(p.sigma));CpR (noshrink_ind,noshrink_ind) = PEB_R(node).Cp;
        q.sigma = CpR;
        KLD_FR(node) = KL(p,q);
    
    end
    
    % Causal Centrality
    CCent = dF_PEB + KLD_FR;
    CCent = CCent (:);
end