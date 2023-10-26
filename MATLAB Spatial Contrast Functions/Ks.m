%spatial contrast
% 2p+1 <= sqrt(N)
% M= area pixel/ area speckle
%Ks needs the function MuInv
function RETURN=Ks(M,N,p)
central=MuInv(M,0,0);

lateral=0;
diagonal=0;
knight=0;
if N>1;
    for eta=1:1:p
        lateral= lateral+(sqrt(N)-eta).*(sqrt(N)).*MuInv(M,eta,0);
        diagonal=diagonal+(sqrt(N)-eta).*(sqrt(N)-eta).*MuInv(M,eta,eta); 
    end
    lateral=4*lateral;
    diagonal=4*diagonal;


    for xi=1:1:p-1;
        for eta=xi+1:1:p;
            knight=knight+(sqrt(N)-eta).*(sqrt(N)-xi).*MuInv(M,eta,xi);
        end
    end
    knight=8*knight;
    correccion=lateral+diagonal+knight;
    correccion=correccion./(N.*(N-1));
else
    correccion=0;
end

RETURN=central-correccion;
RETURN=sqrt(RETURN);
end