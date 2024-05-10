clear all

% Constants (all MKS, except energy which is in eV)
hbar=1.06e-34; q=1.6e-19; qh=q/hbar; a=1e-9;

% inputs
t0=-2.5; D=1e-50; ctr=0; zplus=1i*1e-3;
NL=1; L=zeros(NL); R=L; L(1,1)=1; R(NL,NL)=1;
config=2; %1 for armchair, 2 for zigzag
%NW=floor(14*sqrt(3)); %Armchair
NW=14; %Zigzag

%Hamiltonian
al=t0*[0 1 0 0;1 0 1 0;0 1 0 1;0 0 1 0];
if config==1
    bL=t0*[0 0 0 0;0 0 0 0;0 0 0 0;1 0 0 0];
    bW=t0*[0 0 0 0;1 0 0 0;0 0 0 1;0 0 0 0];
end
if config==2
    bW=t0*[0 0 0 0;0 0 0 0;0 0 0 0;1 0 0 0];
    bL=t0*[0 0 0 0;1 0 0 0;0 0 1 0;0 0 0 0];
end

n=4; %al=4; bW=-1; bL=-1; n=1;
alpha=kron(eye(NW),al)+kron(diag(ones(1,NW-1),+1),bW)+kron(diag(ones(1,NW-1),-1),bW');
%alpha = alpha+kron(diag(ones(1,1),1-NW ),bW ) + kron(diag(ones(1,1),NW-1),bW'); %for CNT's

sigB=zeros(NW*NL*n);siginB=zeros(NW*NL*n);

ii=0; 
for EE=t0*(-0.5:+0.01:+0.5)
ii=ii+1;
ig0=(EE+zplus)*eye(NW*n)-alpha;
if ii==1
gs1=inv(ig0);gs2=inv(ig0);end

BB=0; beta=kron(diag(exp(1i*qh*BB*a*a*(1:1:NW))),bL);
%beta=kron(eye(NW),bL);
H=kron(eye(NL),alpha); if NL>1
H=H+kron(diag(ones(1,NL-1),+1),beta)+kron(diag(ones(1,NL-1),-1),beta');end

change=1;
while change>1e-4
Gs=inv(ig0-beta'*gs1*beta);
change=sum(sum(abs(Gs-gs1)))/(sum(sum(abs(gs1)+abs(Gs))));
gs1=0.5*Gs+0.5*gs1;
end
sig1=beta'*gs1*beta;sig1=kron(L,sig1);gam1=1i*(sig1-sig1');
change=1;
while change>1e-4
Gs=inv(ig0-beta*gs2*beta');
change=sum(sum(Gs-gs2))/(sum(sum(abs(gs2)+abs(Gs))));
gs2=0.5*Gs+0.5*gs2;
end
sig2=beta*gs2*beta';sig2=kron(R,sig2);gam2=1i*(sig2-sig2');

G=inv((EE*eye(NW*NL*n))-H-sig1-sig2);
T(ii)=real(trace(gam1*G*gam2*G'));E(ii)=EE/t0;
if EE==0
T(ii)=T(ii-1); end, EE
end

hold on
h=plot(T,E,'k');
set(h,'linewidth',[1.5])
set(gca,'Fontsize',[10]);
axis([0 10 -0.5 +0.5]);
xlabel('T(E)');
ylabel('(E-\epsilon)/t_0');
title('W = 14*3a')
grid on