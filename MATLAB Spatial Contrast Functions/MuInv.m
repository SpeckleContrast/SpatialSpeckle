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