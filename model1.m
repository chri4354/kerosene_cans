clear all; clc; %close all
set(0,'DefaultFigureWindowStyle','docked')

% for now we're looking at just one potential leader who is outside the group
% (i.e. leader does not form beliefs about other group members;
% and no one else is evaluated as a potential leader)

% expert's belief distribution
aE=4; bE=10; %smaller a, b, larger variance
x=linspace(betainv(.00001,aE,bE),betainv(.99999,aE,bE),100);
funE=@(x)betapdf(x,aE,bE);
E=feval(funE,x);
areaE = integral(funE,min(x),max(x));

% truth belief distribution
aT=10; bT=10;
x=linspace(betainv(.00001,aT,bT),betainv(.99999,aT,bT),100);
funT=@(x)betapdf(x,aT,bT);
T=feval(funT,x);
areaT = integral(funT,min(x),max(x));

% initial individual parameters
aInit=2; bInit=5;

% group size, time, type of network
n=30; tmax=6; net=2;

%%%%%%%%%%%%%%%%%%%%%%%

if net==1; %full net
    neighbas=repmat(1:n,n,1);
elseif net==2; %simple ring lattice, 2 neighbors
    neighbas=[(1:n)+1; (1:n)+2]';
    neighbas(neighbas>n)=neighbas(neighbas>n)-n;
end

figure;
% for each time point
for t=1:tmax
    
    % for each individual
    for i=1:n
        
        % initial individual belief distributions
        if t==1; aI(i,t)=aInit+rand-1; bI(i,t)=bInit+rand-1; end
        
        % current individual belief distribution
        x=linspace(betainv(.00001,aI(i,t),bI(i,t)),betainv(.99999,aI(i,t),bI(i,t)),100);
        funI=@(x)betapdf(x,aI(i,t),bI(i,t));
        I=feval(funI,x);
        areaI(i,t) = integral(funI,min(x),max(x));
        
        %plot
        subplot(3,2,t)
        hold on; plot(E,'r','LineWidth',2)
        hold on; plot(T,'k','LineWidth',2)
        hold on; plot(I,'b')
        
        % overlapIE - overlap between individual and expert's belief distributions
        funIE=@(x)min(feval(funI,x),feval(funE,x));
        overlapIE(i,t)= integral(funIE,0,1);
        
        % overlapET - overlap between true and expert's belief distributions
        funTE=@(x)min(feval(funT,x),feval(funE,x));
        overlapET(i,t)= integral(funTE,0,1);
        
        % evaluation of expert ("respect")
        if t==1; %because need to wait one time point to see evaluations of others
            R(i,t)=mean([overlapIE(i,t) overlapET(t)]);
        else
            R(i,t)=mean([overlapIE(i,t-1) mean(R(neighbas(i,:),t-1)) overlapET(t-1)]);
        end
        
        % adjustment of own belief - parameters adjusted proportional
        % to evaluation of the expert
        aI(i,t+1)=(1-R(i,t))*aI(i,t) + R(i,t)*aE;
        bI(i,t+1)=(1-R(i,t))*bI(i,t) + R(i,t)*bE;
        
    end % end of n
    
    title(strcat('t=',num2str(t)))
end % end of tmax

if net==1; nt='Full'; elseif net==2; nt='Ring'; end
suptitle(strcat(nt,', aE=',num2str(aE),', bE=',num2str(bE), ...
    ', aT=',num2str(aT),', bT=',num2str(bT), ...
    ', aInit=',num2str(aInit),', bInit=',num2str(bInit)))

legend('Expert','True','Individuals')






