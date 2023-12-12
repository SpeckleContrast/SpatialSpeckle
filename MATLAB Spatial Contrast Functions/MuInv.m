% Julio Cesar Juarez-Ramirez, Beatriz Coyotl-Ocelotl, Bernard Choi, Ruben Ramos-Garcia,
% Teresita Spezzia-Mazzocco, Julio C. Ramirez-San-Juan, 
% "Improved spatial speckle contrast model for tissue blood flow imaging: effects of spatial correlation among neighboring camera pixels," 
% J. Biomed. Opt. 28(12) 125002 (5 December 2023) https://doi.org/10.1117/1.JBO.28.12.125002

function RETURN=MuInv(M,eta,xi)
RETURN=(1./((2*pi.*M).^2)).*factorMu(M,eta).*factorMu(M,xi);
end

function RETURN=factorMu(M,indice)
expo=exp(-pi.*M*(indice-1)^2)+...
    -2*exp(-pi.*M*(indice)^2)+...
    exp(-pi.*M*(indice+1)^2);

error=(indice-1)*erf(sqrt(pi.*M)*(indice-1))+...
    -2*(indice)*erf(sqrt(pi.*M)*(indice))+...
    (indice+1)*erf(sqrt(pi.*M)*(indice+1));

RETURN=expo+pi.*sqrt(M).*error;
end
