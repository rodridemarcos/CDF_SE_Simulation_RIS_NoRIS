function R = functionRlocalscattering(M,angleofdeparture,ASDdeg,antennaSpacing)
%Calculates the spatial correlation matrix based on the local scattering
%model.
%
%INPUT:
%M                 = Number of antennas
%angleofdeparture  = Nominal angle of departure (in radians)
%ASDdeg            = Angular standard deviation (in degrees)
%antennaSpacing    = Antenna spacing (in number of wavelengths)
%
%OUTPUT:
%R                 = M x M spatial correlation matrix
%
%
%This Matlab function was developed to generate simulation results to:
%
%Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017),
%"Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency",
%Foundations and Trends in Signal Processing: Vol. 11, No. 3-4,
%pp. 154-655. DOI: 10.1561/2000000093.
%
%For further information, visit: https://www.massivemimobook.com
%
%This is version 1.0 (Last edited: 2017-11-04)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.


%Convert ASD from degrees to radians
ASD = ASDdeg*pi/180;

%Prepare matrix for output
R = zeros(M,M);

%Go through all pairs of antennas
for m1 = 1:M
    for m2 = 1:M

        %Calculate the distance between the antennas
        delta_m = (m1-m2)*antennaSpacing;

        %Calculate the spatial correlation according to (2.15) in the monograph
        R(m1,m2) = besselj(0,2*pi*delta_m*sin(ASD/2));

    end
end

%Rotate the correlation matrix according to the nominal angle of departure
varphi = (0:M-1)'*2*pi*antennaSpacing*sin(angleofdeparture);
rotationVector = exp(1i*varphi);
R = rotationVector*rotationVector'.*R;