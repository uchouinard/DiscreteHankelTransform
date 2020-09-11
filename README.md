# Discrete Hankel Transform

# Matlab code Code for the Discrete Hankel Transform

Previous definitions of a Discrete Hankel Transform (DHT) have focused on methods to approximate the continuous Hankel integral transform without regard for the properties of the DHT itself.  Recently, the theory of a Discrete Hankel Transform was proposed that follows the same path as the Discrete Fourier/Continuous Fourier transform.  This DHT possesses orthogonality properties which lead to invertibility and also possesses the standard set of discrete shift, modulation, multiplication and convolution rules.  The proposed DHT can be used to approximate the continuous forward and inverse Hankel transform.  
 
The full theory can be found in The Discrete Hankel Transform: Properties and Application to the Continuous Hankel Transform, Journal of the Optical Society of America A, Vol. 32, No. 4, pp. 611-622, 2015. http://dx.doi.org/10.1364/JOSAA.32.000611  
 
Description of this code and how to use it can be found in Chouinard U, Baddour N. (2017). Matlab Code for the Discrete Hankel Transform. Journal of Open Research Software. 5(1), p.4. DOI: http://doi.org/10.5334/jors.82 


# Update September 2020
Adi Natan kindly improved some of the code.  
The modifications:
1. The Y matrix code is now vectorized making it ~ x20 faster. 
2. The code has an optional zero padding input similar to Matlab's fft functionality.
3. The code supports also arrays not only vectors similar to Matlab's fft functionality, meaning that the Bessel zeros and Y matrix are only calculated once per run
The new code is posted here as DHT.m and can be found on Adi's github:   https://github.com/adinatan/Discrete-Hankel-Transform
