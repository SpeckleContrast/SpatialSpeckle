% The article "Improved spatial speckle contrast model for tissue blood flow imaging: Effects of spatial correlation among neighboring camera pixels" 
% in the Journal of Biomedical Optics (JBO) should be cited in any work related to the software.

%%%%%%%%%%%%%%%%% 2p+1 <= sqrt(N)  %%%%%%%%%%%%%%
%%%%%%%%    M= area pixel/ area speckle   %%%%%%%
%%%%%%%%   M with exponential separation. %%%%%%%
xmin=1*10^-4;xmax=1*10^1;
Npoints=100;%number of points
M=log(xmin):(log(xmax)-log(xmin))/(Npoints-1):log(xmax);
M=exp(M);
%variables for the plot
TamAx=16;%Axes font size
lw=2;%line width

%Calculating the correlation factors mu, for the first 6 correlation
%factors
Mu00=MuInv(M,0,0);
Mu10=MuInv(M,1,0); 
Mu11=MuInv(M,1,1);
Mu20=MuInv(M,2,0);
Mu21=MuInv(M,2,1);
Mu22=MuInv(M,2,2);

figure('Renderer','painters','Position',[30,60,1400,500]);
subplot(1,2,1);
semilogx(M,Mu00,'m','LineWidth',lw);hold on;
semilogx(M,Mu10,'r','LineWidth',lw);
semilogx(M,Mu11,'k','LineWidth',lw);
semilogx(M,Mu20,'g','LineWidth',lw);
semilogx(M,Mu21,'b','LineWidth',lw);
semilogx(M,Mu22,'c','LineWidth',lw);hold off;
grid;
legend('1/\mu_{0,0}','1/\mu_{1,0}','1/\mu_{1,1}',...
    '1/\mu_{2,0}','1/\mu_{2,1}','1/\mu_{2,2}','Location','southwest');
xlabel('M = a^2/\pi b^2 = Pixel area/Speckle area'); ylabel('1/\mu_{\eta,\xi} (a.u.)');
xlim([xmin, xmax]);
set(gca,'fontsize',TamAx);


%contrast for 1x1, 3x3, 5x5 sliding window respectively
K11=Ks(M,1^2,0);
K33=Ks(M,3^2,1);
K55=Ks(M,5^2,2);

subplot(1,2,2);
semilogx(M,K11,'b','LineWidth',lw);hold on;
semilogx(M,K33,'m','LineWidth',lw);
semilogx(M,K55,'k','LineWidth',lw);
hold off;
grid;

legend('K_{s}(1,0) (current model)','K_{s}(3^2,1)','K_{s}(5^2,2)','Location','southwest');
xlabel('M = a^2/\pi b^2 = Pixel area/Speckle area'); ylabel('K_{s}(N,p) (a.u.)');
ylim([0.3,1]);xlim([0.01,10]);
set(gca,'fontsize',TamAx);