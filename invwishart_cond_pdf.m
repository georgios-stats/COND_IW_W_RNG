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

% Compute the density of the inverse Wishart distribution 
% conditional on the 1st block diagonal sub-matrix.

% Karagiannis, G., Andrieu, C. (under revision). Drawing samples from 
% inverse Wishart distributions conditioning on the 1st block diagonal 
% sub-matrix; with an application to variable selection on a GLMM model 
% with nested random effects

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
%     X : The complete (p_1+p_2)X(p_1+p_2) semi-positive matrix.
%     p_1: The 1st dimention of the 1st block diagonal sub-matrix.
%     p_2: p_2 = size(X,1)-p_1 .
%     m : degrees of freedom of the big matrix
%     S : The complete (p_1+p_2)X(p_1+p_2) semi-positive scale matrix.
% Outputs :
%     pdf: The value of density
%     logpdf:  log(pdf)
    
function [pdf,logpdf] = invwishart_cond_pdf(X, p_1, p_2, m, S)

    [~, AA] = invwishpdf(X, p_1+p_2, m, S) ; 
    [~, BB] = invwishpdf(X(1:p_1,1:p_1), p_1, m-p_2, S(1:p_1,1:p_1)) ; 

    logpdf = AA -BB ;
    pdf = exp(logpdf) ;

end


