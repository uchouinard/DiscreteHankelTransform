%{

    This script demonstrates how to use the discrete hankel transform 
    

    The discrete transform is performed on a set of N-1 points 

    The hankel transform is of order "n"

    Two cases are demonstrated in this script: the forward and inverse
    transform

_______________________________________________________________________
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
________________________________________________________________________
%}

%the order and number of sample points will be taken as the same for both
%cases, the bessel zeros will also be used in both cases

n=1; %order
N=256; %number of points

%the first step is to create a vector of bessel zeroes that will be used
%at multiples places in the code with "besselzero"

%the last entry in besselzero is the kind of the bessel function, 1 is used
zeros=besselzero(n,N,1);

%both forward and inverse transformed are performed with the same
%Transformation matrix "yMatrix"

yMatrix=YmatrixAssembly(n,N,zeros);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Case 1: Forward Transform
%}


% let's define the space limitation. 
R=15; 


%the function that has to be transformed needs to be discretized at
%specific sample points. this can be achieved by using the spaceSampler
%function

samplePointsSpaceDomain=spaceSampler(R, zeros);


% once the sample points are defined, the function is discretized
%the function that will be used is the sinc function defined as:
% f(r)=sin(a*r)/(a*r)

%let's take an arbitrary factor a;
a=5;


%then the discrete function is:
f=sin(a*samplePointsSpaceDomain(:))./(a*samplePointsSpaceDomain(:));


%in order to perform the transform, the discrete function has to be a 
%column vector 

%there is also the need to apply a scaling factor defined as:

scalingFactor=R^2/zeros(N); %space limited function scaling factor


%the transform is then performed by:
Ftransformed=yMatrix*f(:)*scalingFactor;

%the corresponding frequency domain points are given by:

samplePointsFrequencyDomain=freqSampler(R,zeros);

%and the result can be seen as:

figure()
plot(samplePointsFrequencyDomain,Ftransformed);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
 Case 2: Inverse Transform
%}


%in the case of the inverse transform, the function is limited in the
%frequency domain by "W"

W=30; 

%the sampling is then performed using the "freqSampler" function. However
%the frequency limit "W" has to be transformed to "R" by:

R=zeros(N)/W;

frequencySamplePoints=freqSampler(R,zeros);


%the discrete function is then given by:

F=  (frequencySamplePoints<a).*(((frequencySamplePoints./a).^n).*cos(n*pi/2)./ ...
    (a^2.*sqrt(1-frequencySamplePoints.^2./a^2).*(1+sqrt(1-frequencySamplePoints.^2/a^2)).^n))+ ...
    (frequencySamplePoints>a).*((sin(n.*asin(a./frequencySamplePoints)))./(a^2.*sqrt(frequencySamplePoints.^2./a^2-1)));
    
%the scaling factor for a frequency limited function while performing a
%IDHT is given by:


scalingFactorFreq=1;

%the inverse transform is then given by:

fTransformed=yMatrix*F(:)*scalingFactorFreq;


%and the corresponding sample points in the space domain are given by:

spaceSamplePoints=spaceSampler(R,zeros);


%the result can be seen by:

figure()
plot(spaceSamplePoints,fTransformed)

%% Compare DHT computed results with continuous function
%Check results for forward transform
figure()
plot(spaceSamplePoints,fTransformed,samplePointsSpaceDomain,f)
legend('IDHT computed function','continuous function')

%Check results for inverse transform
figure()
plot(samplePointsFrequencyDomain,Ftransformed,frequencySamplePoints,F);
legend('DHT computed function','continuous function')
