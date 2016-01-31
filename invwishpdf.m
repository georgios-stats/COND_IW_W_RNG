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
% with expected value E(X) = S/(m-p-1)

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


function [ pdf, logpdf ] = invwishpdf( X, p, m, S ) 

    if ( m <= p-1.0 )
        fprintf('ERROR::WISHPDF::(m <= p-1.0)') ;
    end

    logMVG_value = 0.25*p*(p-1.0)*log(pi) ...
                    +sum(gammaln(0.5*(m+1-(1:p)))) ;

    logpdf = -0.5*m*p*log(2.0) ...
                +0.5*m*log(det(S)) ...
                -logMVG_value ...
                +0.5*(m+p+1)*log(det(X)) ...
                -0.5*trace(S/X) ;

    pdf = exp( logpdf ) ;
    
end




