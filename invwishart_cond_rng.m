% Copyrigtht 2012, 2016 Georgios Karagiannis
%
% This file is part of COND_IW_W_RNG.
%
% COND_IW_W_RNG is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
%
% COND_IW_W_RNG is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with COND_IW_W_RNG.  If not, see <http://www.gnu.org/licenses/>.

% Generate random values from the inverse Wishart distribution 
% conditional on the 1st block diagonal sub-matrix.

% Karagiannis, G., Andrieu, C. (2009, 2016). Drawing samples from 
% inverse Wishart distributions conditioning on the 1st block diagonal 
% sub-matrix; with an application to variable selection on a GLMM model with 
% nested random effects, GitHub repository manuscript

% Georgios Karagiannis
% 
% Postdoctoral research associate
% Department of Mathematics, Purdue University
% 150 N. University Street
% West Lafayette, IN 47907-2067, USA
% 
% Telephone: +1 765 496-1007
% 
% Email: gkaragia@purdue.edu
% 
% Contact email: georgios.stats@gmail.com


% Inputs :
%     p_1: The 1st dimention of the 1st block diagonal sub-matrix.
%     p_2: p_2 = size(X,1)-p_1 .
%     m : degrees of freedom of the big matrix
%     Sig_11 : The p_1Xp_1 1st block diagonal sub-matrix we condition.
%     S : The complete (p_1+p_2)X(p_1+p_2) semi-positive scale matrix.
% Outputs :
%     Sig: The complete random (p_1+p_2)X(p_1+p_2) semi-positive matrix.
%     cholSig:  chol(Sig,'lower')



function [Sig, cholSig] = invwishart_cond_rng(m, p_1, p_2, Sig_11, S) 

    C_11 = chol(Sig_11,'lower') ;
    Sig = chol(S,'lower') ;
    L_11 = Sig(1:p_1,1:p_1) ; 
    L_21 = Sig((p_1+1):(p_1+p_2),1:p_1) ;
    L_22 = Sig((p_1+1):(p_1+p_2),(p_1+1):(p_1+p_2)) ;

    N_21 = randn(p_2,p_1) ;
    
    B_22 = zeros(p_2,p_2) ;
    for i = 1:p_2
        for j = 1:p_2
            if (j>i)
                B_22(i,j) = randn() ;
            elseif (j==i)
                %B_22(i,j) = sqrt( sum( ( randn((m-p_2+i),1) ).^2 ) ) ;
                %B_22(i,j) = sqrt( chi2rnd(m-p_2+i) ) ;
                B_22(i,j) = sqrt( 0.5*randg(0.5*(m-p_2+i)) ) ;
            end
        end
    end

    B_22 = L_22/(B_22') ;
    cholSig = [C_11 , zeros(p_1,p_2) ; (L_21-B_22*N_21)*(L_11\C_11) , B_22] ;
    
    Sig = cholSig*(cholSig') ;

end

