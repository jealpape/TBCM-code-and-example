clc; close all; clear;
%% Example of how to use the TBCM 
%Name: José Manuel Rojas // Jes�s Parra
%Date: 19/07/22
%addpath("..\..\TBCM\")


addpath("tools/")
addpath("PhonationModelsCode/")

resul = './Prueba/';
vowels = {'AE', 'A'};
gender = {'male', 'female'};
fs = 44100;
ig = 1;

Pr = 1100;
a_PCA = 0; %0:0.1:0.1;
ACT = combvec(Pr,a_PCA);

vPL = ACT(1,:);
vPCA = ACT(2,:);
D = length(vowels)*length(vPL)*length(gender);

a_CT = 0.2; %0:0.1:1; 
a_TA = 0.4; %0:0.1:1; 
a_LCA = 0.6; %0.2:0.1:0.8;
ACT = combvec(a_CT,a_TA,a_LCA);

vCT = ACT(1,:);
vTA = ACT(2,:);
vLCA = ACT(3,:);
N = length(vCT);

clear a_CT a_TA a_LCA a_PCA Pr ACT 

%%   

%poolobj = parpool('local',2); % length(vowels)
  
for n = 1:N 
    
    SimResults = struct();
    
    ct = vCT(n);
    ta = vTA(n);
    lca = vLCA(n);
    ia = lca;

    H1H2 = zeros(D,1); HRF = zeros(D,1); MFDR = zeros(D,1);
    ACFL = zeros(D,1); SQ = zeros(D,1); OQ = zeros(D,1); CPP = zeros(D,1);
    F0 = zeros(D,1); SPL = zeros(D,1); PS = zeros(D,1); OQN = zeros(D,1);
    PC = zeros(D,1); PL = zeros(D,1); 
    
    CT = zeros(D,1); TA = zeros(D,1); LCA = zeros(D,1);  PCA = zeros(D,1);
    PGO = zeros(D,1); VOW = zeros(D,1); GEN = zeros(D,1);
    
    id = 1;
    for a = 1:length(vowels)                %PCA activation
        for b = 1:length(vPL)
            for c = 1:length(gender)
                
                pl = vPL(b);
                pca = vPCA(b);
                
                act = [lca;lca;0;ct;ta];

                [dX,Lg,X,a_g,Ut,Pgo,Pout,Psub,Pcol,tm] = SimuBaseSignal(act,pl,vowels{a},gender{c});

                [asa, asp, eje, h1h2, hrf, mfdr, acfl, f0, sq, oq, oqn, cpp, naq] = features(X,X,Ut,fs);
                spl=measures_getSPL(Pout,fs); ps=mean(Psub); pgo = Pgo(end);
                [P_col_max, ~]=findpeaks(Pcol);
                pc=mean(P_col_max);
                
                
                % save signals
                SimResults.dX(:,:,id)=dX; SimResults.Lg(:,id)=Lg; SimResults.X(:,:,id)=X;
                SimResults.a_g(:,id)=a_g; SimResults.Ut(:,id)=Ut; SimResults.PGO(:,id)=Pgo;
                SimResults.Pout(:,id)=Pout; SimResults.Psub(:,id)=Psub; SimResults.Pcol(:,id)=Pcol; SimResults.act(:,id)=act;
                SimResults.PL(id)=pl;  SimResults.gender{id}=gender{c};  SimResults.vowel{id}=vowels{a};

                % save table
                GEN(id)= c; VOW(id) = a; PL(id) = pl; CT(id) = ct; TA(id) = ta; LCA(id) = lca; PCA(id) = pca; 
                PGO(id) = pgo; H1H2(id) = h1h2; HRF(id) = hrf; MFDR(id) = mfdr; ACFL(id) = acfl;
                SQ(id) = sq; OQ(id) = oq; OQN(id) = oqn; F0(id) = f0; SPL(id) = spl; PS(id) = ps; PC(id) = pc; CPP(id) = cpp;

                % disp(['simu ' int2str(id)]);
                id = id+1;
            end
        end
    end
    
    vec = [GEN, VOW, PCA, CT, TA, LCA, PL, H1H2, HRF , MFDR, ACFL, ...
           SQ, OQ , OQN, F0, CPP, SPL, PS, PC, PGO];
    
    Tab= array2table (vec, 'VariableNames', {'GEN' 'VOW' 'PCA' 'a_CT' 'a_TA' 'a_LCA' 'PL' 'H1H2' 'HRF' 'MFDR' 'ACFL' 'SQ' 'OQ' 'OQN' 'F0' 'CPP' 'SPL' 'PS' 'PC' 'PGO'});
    
    parsave([resul 'Table/Tab_' int2str(n)],Tab)
    parsave([resul 'Signal/Data_' int2str(n)],SimResults)
    
    
    
end
%%
%delete(poolobj)
