%%
% CollisionForces: Function implementing the computation of the 
% collision forces produced during closure phase in accordance to the  
% triangular body-cover model of the vocal folds.
%
% Structure: [Fu_col,Fl_col] = CollisionForces(TBCMobj)
% where
%
% BCMObj: is an object from BodyCoverModel (handle) class,
% Fu_col: is the collision force in the upper mass,
% Fl_col: is the collision force in the lower mass.
%
% References:
% [1] G. E. Galindo, S. D. Peterson, B. D. Erath, C. Castro, R. E. Hillman,
%     and M. Zañartu, “Modeling the Pathophysiology of Phonotraumatic Vocal
%     Hyperfunction With a Triangular Glottal Model of the Vocal Folds,” 
%     J. Speech Lang. Hear. Res., vol. 60, no. 9, p. 2452, Sep. 2017.
% [2] P. Birkholz, B. J. Kröger, and C. Neuschaefer-Rube, “Articulatory 
%     synthesis ofwords in six voice qualities using a modified two-mass 
%     model ofthe vocal folds.” Paper presented at the First International
%     Workshop on Performative Speech and Singing Synthesis, Vancouver,
%     British Columbia, Canada, 2011.
% [3] P. Birkholz, B. J. Kröger, and C. Neuschaefer-Rube, “Synthesis of 
%     breathy, normal, and pressed phonation using a two-mass model with a 
%     triangular glottis.” In P. Cosi, R. De Mori, G. Di Fabbrizio, & R.
%     Pieraccini (Eds.), Interspeech 2011: 12th Annual Conference of the
%     International Speech Communication Association (pp. 2681–2684).
%     Baixas, France: Interna- tional Speech Communication Association. 2011   
%
% Coded by Gabriel Alzamendi, January 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.

% IMPORTAN!! 
% There are two typos in [1] in the Eqs. (A19) and (A20)
% (A19) ... + ñ_{uCol} (x_u^3 + 1.5 x_u^2 \deltax_u \alpha_u + x_u \deltax_u^2 [+] \alpha_u^2 + 0.25 \deltax_u^3 \alpha_u^3)    -->>(CHANGE FOR)-->>
%       ... + ñ_{uCol} (x_u^3 + 1.5 x_u^2 \deltax_u \alpha_u + x_u \deltax_u^2 \alpha_u^2 + 0.25 \deltax_u^3 \alpha_u^3)  
% (A20) ... + ñ_{lCol} (x_l^3 + 1.5 x_l^2 \deltax_l \alpha_l + x_l \deltax_l^2 [+] \alpha_l^2 + 0.25 \deltax_l^3 \alpha_l^3)    -->>(CHANGE FOR)-->>
%       ... + ñ_{lCol} (x_l^3 + 1.5 x_l^2 \deltax_l \alpha_l + x_l \deltax_l^2 \alpha_l^2 + 0.25 \deltax_l^3 \alpha_l^3)  


function [Fu_col,Fl_col] = CollisionForces(TBCMobj)
  % Mass displacements
  xu = TBCMobj.xData(1);
  xl = TBCMobj.xData(2);
  xu_col = TBCMobj.xu_col;
  xl_col = TBCMobj.xl_col;
  deltau0 = 0.5*max([0, TBCMobj.xi_02]);  % Posterior position of the upper mass in the TBCM. 
  deltal0 = 0.5*max([TBCMobj.xi_01-TBCMobj.xi_02, TBCMobj.xi_01]); % Posterior position of the lower mass in the TBCM.
  alphau = TBCMobj.alpha_u;
  alphal = TBCMobj.alpha_l;
  
  xiu_f = TBCMobj.xData(1)-deltau0-TBCMobj.xu_col;
  xil_f = TBCMobj.xData(2)-deltal0-TBCMobj.xl_col;
  
  x_refu = xiu_f + alphau*abs(deltau0);
  x_refl = xil_f + alphal*abs(deltal0);
  
  % Collision forces
  Fu_col = - TBCMobj.hu_col*alphau*( x_refu ...                 % [N] Collision force in the upper mass   *TBCMobj.Lg
              + TBCMobj.NonLinMode*TBCMobj.etau_col*x_refu*(x_refu^2+(alphau*deltau0)^2)  );

  Fl_col = - TBCMobj.hl_col*alphal*( x_refl ...                 % [N] Collision force in the lower mass    *TBCMobj.Lg
              + TBCMobj.NonLinMode*TBCMobj.etal_col*x_refl*(x_refl^2+(alphau*deltau0)^2)  );  
end

