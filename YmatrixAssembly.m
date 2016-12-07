% Y is the N-1 x N-1 transformation matrix to be assembled
%
% n is the order of the bessel function
% N-1 is the size of the transformation matrix
%

%{


Author : Ugo Chouinard
University of Ottawa
March 2015

Copyright (c) 2015, Ugo Chouinard
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

%}


function Y = YmatrixAssembly(n,N,zeros)


for m=1:N-1
    
    for k=1:N-1
        
        % the zeros that are used to compute the matrix
        %they are retrieved from the zero array for clarity
        jnk=zeros(k);
        jnm=zeros(m);
        jnN=zeros(N);
       
        %bessel function of order n+1
        jnplus1=besselj(n+1, jnk);
        
        Y(m,k)=(2*besselj(n,(jnk*jnm/jnN)))/(jnN*jnplus1^2);
        
        
    end
end



end

