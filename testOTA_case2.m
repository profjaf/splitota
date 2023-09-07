%this code applies OTA for simulated data or data obtained from 
%simulation of constant-amp orders
%three orders: 1X, 1.2X and 40 Hz
% Case-1: 5 to 60 Hz, duration: 30 sec, fs = 1280, bnr= rev
clear;
close all;
addpath("functions\");
caseNo = 2; %1: constant-amp, 2: two-DOF
thetSrc = 0; %angular position source, 0: actual, 1: corrected by interpolation of speed, 2: uncorrected; 
resSrc = 0; %response source, 0: simulation, 1: SSSS plate matlab, 2: plate Ansys
opt = 1; %1:vkf, 2:svkf

%rmmat=[2 5 10 20 40 60 80 100 120 140 160 180 200 250 300 350 400 450 500]; %VKF weighting factor, 
rmmat=[43];
rmSplitmat=sqrt(0.5)*rmmat;
forder=0;  %VK filter order: 0,1,2, No. of poles = filter order + 1

fs = 1280;
dwnsamp = 1;%downsampling factor of vibration data, 1: no downsampling
bsize = 200;%processing length for constant DT in samples, 4800
overlap = 0;%overlapping for constant DT
bnr = 4; %processing block in rev
nroverlp=0;

%define orders; 
% maximum order frequency MUST NOT exceed half sampling rate
% Orders must be in ascending order for Run-up test and descending order
% for coast-down test
% constant frequency items can be put in any location, better at the
% end of orders
% OrderType: 0: Order, 1: constant frequency
NTraces = 3;
TraceOrder= zeros(NTraces,1);
TraceType=zeros(NTraces,1);

TraceOrder(1)=1;
TraceType(1)=0;

TraceOrder(2)=2;
TraceType(2)=0;

TraceOrder(3)=40;
TraceType(3)=1;

TrigLev1=0;
HystV = 0.2;

TotalTime = 30;
Np = TotalTime*fs;
Nst=800;Nend = Np-200; %trim edges
Res2 = zeros(Np,2);
ResCF = zeros(Np,1);
ft = zeros(Np,1);


dt = 1/fs;
ts=dt*(0:Np-1)';
freqHz=5;
thetaM = zeros(Np,1);
f1 = zeros(Np,1);
f3 = TraceOrder(3)*ones(Np,1);
ang=0;
%% define the simulation for two-DOF system
K = 3.0*[30000 -16000;-16000 26000];
M =[1.2 0; 0 1.0];
C=[18 0; 0 15];
[phi,D]=eig(K,M); % modes shapes and eigenvalues, 
%               Note: matlab provides mass normalized mode shapes
wn=sqrt(D);
z1=(C(1,1)/M(1,1))/(2*wn(1,1));
z2=(C(2,2)/M(2,2))/(2*wn(2,2));
wr=[wn(1,1) wn(2,2)];
zr=[z1 z2];
alphar=@(i,j,r,w) 10000*phi(i,r)*phi(j,r)/(wr(r)*wr(r)-w*w+2i*zr(r)*w*wr(r));
%% start signal generation
acta = zeros(Np,3);
actb = zeros(Np,3);
actamp = zeros(Np,3);
% f1=5+(55/(Np-1))*(0:Np-1)';
for i=1:Np
%     ang = 2*pi*(5.0*ts(i)+55.0*ts(i)*ts(i)/(2*TotalTime));%direct integration of omega
    w = freqHz*2*pi;
    if caseNo == 1
       for j=1:3
           if j == 1
               actb(i,j) = 0;%4 + 4* ((i-1)/(Np-1))^2;
               acta(i,j) = 4 + 4* ((i-1)/(Np-1))^2;
           elseif j==2
               actb(i,j) = 0;%2 + 2* (i-1)/(Np-1);
               acta(i,j) = 2 + 2* (i-1)/(Np-1);
           else
               actb(i,j) = 0;%3;
               acta(i,j) = 3.5;
           end
       end
    else
      a11_1X=alphar(1,1,1,w)+alphar(1,1,2,w);
      acta(i,1)=real(a11_1X);
      actb(i,1)=imag(a11_1X); 
      a11_4X=alphar(1,1,1,TraceOrder(2)*w)+alphar(1,1,2,TraceOrder(2)*w);
      acta(i,2)=real(a11_4X);
      actb(i,2)=imag(a11_4X);       
      acta(i,3)=sqrt(0.2);
      actb(i,3)=sqrt(0.2); 
    end
    for j=1:3
         actamp(i,j) = sqrt(acta(i,j)^2+actb(i,j)^2);
    end
    Res2(i,1) = acta(i,1)*cos(TraceOrder(1)*ang) + actb(i,1)*sin(TraceOrder(1)*ang) ...
        + acta(i,2)*cos(TraceOrder(2)*ang) + actb(i,2)*sin(TraceOrder(2)*ang) ...
        + acta(i,3)*cos(2*pi*TraceOrder(3)*ts(i)) + actb(i,3)*sin(2*pi*TraceOrder(3)*ts(i));
    
    f1(i) = freqHz;    
    Res2(i,2) = 5*sin(ang);
    thetaM(i) = ang;
    ang = ang + w*dt;%digital integration of Omega
    freqHz = 5 + 65 *i/(Np-1); % = 5 + 55 t/T
end
% thetaM=2*pi*cumsum(f1,1)*dt;
% plot(env1x)
% return
f2 = f1*TraceOrder(2);
r = Res2(:,2); %last column is tacho signal

if(resSrc == 1)
   filePath = 'G:/My Documents/MATLAB/ModalAna/pLSCF/input/';
   Res3 = load([filePath,'timeDataPlateOTA2.txt']);
   Res2(:,1) = 1000*Res3(:,1);%read from response file, most point affected by uncorrected theta is 12
   r = Res3(:,13);
end

if(resSrc == 2)
   filePath = 'E:/plate_CFFF/';
   Res3 = load([filePath,'acc1noCF.txt']);
   r1 = rms(Res3);%calculate rms level of the signal in order to add specified noise
   noiselev = 0.00*r1;% gives 20log(1/0.125)
   Res3 = Res3 +noiselev*randn(length(Res3),1);
   Res2(:,1) = Res3/1000 + 1*ResCF(:,1); %repelem(Res3,3);%upsample since Ansys use sample rate of 2000 instead of 6000
                              %due to large result file

end
% insig=Res2(1:360,1);
% plot(insig)
% return

L = size(Res2,1);
firstdetected = 0;
theta = zeros(L,1);%extracted from tacho signal
thetaJ = zeros(L,1);%extracted from poly fitting
cpos = zeros(L,1);
OMG = zeros(L,1);
OMG0 = zeros(L,1); %from unfiltered tacho signal
pindex0=1;

nr=0;
Jterm=0;
stp=0; %starting or first edge index
endp=L; %last edge index
if r(1) > TrigLev1 
    PState=1;
else
    PState=0;
end
for i=2:L
    if r(i) > (TrigLev1 + HystV)
        if PState == 0
            if firstdetected > 0
                Count1 = i - pindex0 ;
                endp = i;
                nr = nr+1; %number of revoultions
                cpos(nr) = i;

                t2=  Count1;
                if nr > 1
                    w1= 2*pi*fs / t1; %assuming linear speed variation, calculate the slope
                    w2=2*pi*fs / t2;
                    a1 = (w2-w1)/t2;%dw/dt
                else
                    w1 = 2 * pi *fs / Count1;
                    w2 = w1;
                    a1 = 0;
                end
                tt=0;
                for j = 1:Count1
                    w = w1 + a1*j; %instantaneous speed
                    tt = tt + w/fs; %to ensure complete 2pi per cycle
                end
                Jterm = w1 - 2 * pi *fs / (Count1 + 2*randn(1,1));
                for j = 1:Count1
                    w = w1 + a1*j; %instantaneous speed
                   % OMG0(j+pindex0) = (w1+w2)/2; %good estimate
                    OMG0(j+pindex0) = w1+(w2-w1)*j/Count1; %good estimate
                    theta(j+pindex0) = theta(j-1+pindex0) + (w*dt)*2*pi/tt;
                   
                   % theta(j+pindex0) = theta(j-1+pindex0) + (w2*dt)+0.000002*i/L;
                end
                psi0 = theta(Count1+pindex0);
                
                t1= Count1;
                pindex0 = i;
            else
                stp = i;
                firstdetected = 1;
                pindex0 = i;
                theta(pindex0) = 0.05;
               
                nr = 0;                
                t1 = 0 ;
            end
        end
        PState = 1;
    end
    if r(i) < (TrigLev1 - HystV)
        PState = 0;
    end
end
rf=zeros(1,L);
rev=2;
for i = 1:L
    if i <= cpos(1)
        r(i)= 0;%r(i)/5;
    end
    if cpos(rev) == i %update filter parameters each rev
        i0 = cpos(rev-1);       
       
        ci=cpos(rev)-i0;
        for j=1:ci
            r(i0+j-1)=sin((j-1)*2*pi/ci); %this gives better results
        end
        rev = rev+1;
    end
end


fncr= zeros(nr+12,1);
fncrm= zeros(nr+12,1);
cposm= zeros(nr+12,1);
for i=2:nr
    fncr(i) = fs / (cpos(i) - cpos(i-1));%just to show the accuracy of the speed signal extracted      
end
for i=2:nr-1
    fncrm(i) = (fncr(i) + fncr(i+1))/2;       
end
for i=2:nr   
    cposm(i)=(cpos(i)+cpos(i-1))/2;
end
fncrm(nr)=fncr(nr);
figure(1)
%plot(cpos(2:rev),fncr(2:rev),'r');
hold on
%plot(cpos2(2:rev),fcr(2:rev),'b')
hold on
%yy = csaps(cpos(2:rev),fncr(2:rev),0,cpos(2:rev));
fspd=fit(dt*cpos(2:nr),fncrm(2:nr),'poly3','Robust','Bisquare') ;%poly3, try smoothingspline ,'Robust','Bisquare' ,'LAR'
%fspd=fit(dt*cposm(2:nr),fncr(2:nr),'poly3','Robust','Bisquare') ;%poly3, try smoothingspline ,'Robust','Bisquare' ,'LAR'
yy=fspd(dt*(0:L-1));
plot(dt*(0:L-1),yy)
hold on
plot(dt*(0:L-1),f1);
thetaJ=2*pi*cumsum(yy,1)*dt;
% thetaJ=thetaJ-thetaJ(1);
% thetaJ=[0;thetaJ(1:Np-1)]; %insert zero at the begining since, angle is zero at startup



fprintf('the number of revolutions:  %d \n',nr);
%  return
figure(2)
theterr=thetaM-thetaJ;
plot(theterr(1:Np-1000))
    %return
if(thetSrc == 0)
    theta2=thetaM;%use actual
end
if(thetSrc == 1)
    theta2=thetaJ;%use interpolated
end
if(thetSrc == 2)
    theta2=theta;% use un-corrected theta
    %theta2= theta2+0.03*randn(Np,1);
end
%   return
%% begin OTA
No = 1;
rn=length(rmmat);
mpmat=zeros(rn,3,2);
for opt=1:1
    for kk=1:rn
        if opt == 1 % vkf
            %      Np=2000;
            clf
            rm=[rmmat(kk) rmmat(kk) rmmat(kk)];
            FRFOrder = zeros(Np,NTraces);

            blktime = dt*(0:Np-1);

            theta3=2*pi*cumsum(f3,1)*dt;
            theta3=theta3-theta3(1);
            if(thetSrc == 0)
                vktheta=[thetaM TraceOrder(2)*thetaM theta3];%use actual
            end
            if(thetSrc == 1)
                vktheta=[thetaJ TraceOrder(2)*thetaJ theta3]; % use corrected frequencies
            end
            if(thetSrc == 2)
                vktheta=[OMG0 TraceOrder(2)*OMG0 theta3]; % use un-corrected frequencies
            end
            xe=Res2(1:dwnsamp:Np,1);
            [xm,bwm] = vkmmy(xe,vktheta,fs,rm,forder);

            FRFOrder(:,1) = xm(:,1);
            FRFOrder(:,2) = xm(:,2);
            FRFOrder(:,3) = xm(:,3);

            figure(3)
            plot(blktime(Nst:Nend),abs(FRFOrder(Nst:Nend,1)));%ylim([0 max(abs(FRFOrder(Nst:Nend,1,1)))])
            hold on
            plot(blktime(Nst:Nend),actamp(Nst:Nend,1));ylim([0 max(abs(actamp(Nst:Nend,1)))])
            figure(4)
            plot(blktime(Nst:Nend),abs(FRFOrder(Nst:Nend,2)));ylim([0 max(abs(FRFOrder(Nst:Nend,2)))])
            figure(5)
            plot(blktime(Nst:Nend),abs(FRFOrder(Nst:Nend,3)));ylim([0 max(abs(FRFOrder(Nst:Nend,3)))]);

            mytheta = exp(1i*vktheta);
            reconssig = zeros(Np,1);
            for i=1:Np
                reconssig(i) = xm(i,1)*mytheta(i,1) + xm(i,2)*mytheta(i,2) + xm(i,3)*mytheta(i,3);
            end
            %     figure(6)
            %     plot(imag(reconssig));
            %     title('Imag. of Reconstructed signal', 'FontSize', 11, 'Color', 'k')
            %     ylabel({'Imag.'}, 'FontSize', 10);
            %     figure(7)
            %     plot(real(reconssig));
            %     title('Real of Reconstructed signal', 'FontSize', 11, 'Color', 'k')
            %     ylabel({'Real'}, 'FontSize', 10);
            %     figure(8)
            %     plot(xe-real(reconssig));
            %     title('Differnce actual - real(Reconstructed)', 'FontSize', 11, 'Color', 'k')
            %     ylabel({'Diff.'}, 'FontSize', 10);

            figure(9)
            bsize=2;
            endplot=fix(Nend/bsize);
            bp=fix(2000/bsize);
            plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,1)),'-bs','LineWidth',1,'MarkerEdgeColor','b', 'MarkerSize',4,'MarkerIndices',1:bp:endplot);ylim([0 max(1.1*abs(FRFOrder(Nst:bsize:Nend,1,1)))]);
            hold on
            plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,2)),'-rd','LineWidth',1,'MarkerEdgeColor','r', 'MarkerSize',4,'MarkerIndices',1:bp:endplot)
            plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,3)),'-go','LineWidth',1,'MarkerEdgeColor','g', 'MarkerSize',4,'MarkerIndices',1:bp:endplot)
            plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,1)),':b','LineWidth',1.2);ylim([0 1.1*max((actamp(Nst:Nend,1)))]);
            plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,2)),':r','LineWidth',1.2);
            plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,3)),':g','LineWidth',1.2);
            set(gcf,'position',[400,400,400,300])
            legend('1X',append(num2str(TraceOrder(2)),'X'),append(num2str(TraceOrder(3)),'Hz'),'actual');
            ylabel('Amp. (mm)');xlabel('time (s)')

            MAEa=sum(abs(acta(Nst:Nend,:)-real(FRFOrder(Nst:Nend,:))),1)/(Nend-Nst)
            MAEb=sum(abs(actb(Nst:Nend,:)+imag(FRFOrder(Nst:Nend,:))),1)/(Nend-Nst)
            MPEA=100*max(abs(actamp(Nst:Nend,:)-abs(FRFOrder(Nst:Nend,:))))./max(actamp(Nst:Nend,:))
            mpmat(kk,:,opt)=MPEA;
            %MPEA=100*max(abs(actamp(Nst:Nend,:)-abs(FRFOrder(Nst:Nend,:)))./actamp(Nst:Nend,:))

        end
        if opt == 2 % vk split
            clf
            rmSplit=[rmSplitmat(kk) rmSplitmat(kk) rmSplitmat(kk)];
            FRFOrder = zeros(Np,NTraces);
            xa=zeros(Np,NTraces);
            xb=zeros(Np,NTraces);
            blktime = dt*(0:Np-1);
            theta3=2*pi*cumsum(f3,1)*dt;
            theta3=theta3-theta3(1);
            if(thetSrc == 0)
                vktheta=[thetaM TraceOrder(2)*thetaM theta3];%use actual
            end
            if(thetSrc == 1)
                vktheta=[thetaJ TraceOrder(2)*thetaJ theta3]; % use corrected frequencies
            end
            if(thetSrc == 2)
                vktheta=[OMG0 TraceOrder(2)*OMG0 theta3]; % use un-corrected frequencies
            end

            shift=0;
            xe=Res2(1+shift:dwnsamp:Np+shift,1);
            % [xm,bwm] = vkmmy(xe,fp,fs,rm,forder); this and vkNeilOneOrd provide
            % same results for single order tracking
            [xm,bwm] = vkSplitMultiOrd(xe,vktheta,fs,rmSplit,forder);
            % [xm,bwm] = vkDirctOneOrd(xe,fp,fs,0.01,forder);
            for i=1:NTraces
                xa(:,i)=xm(1:Np,i);
                xb(:,i)=xm(1+Np:end,i);
                FRFOrder(:,i) =sqrt(xa(:,i).*xa(:,i)+xb(:,i).*xb(:,i));
            end
            for selOrd=1:3
                figure(2+selOrd)
                plot((FRFOrder(Nst:Nend,selOrd)));ylim([0 1.1*max(abs(FRFOrder(Nst:Nend,selOrd)))])
                hold on
                plot(actamp(Nst:Nend,selOrd));ylim([0 1.1*max(abs(actamp(Nst:Nend,selOrd)))])
            end

%             reconssig = zeros(Np,1);
%             for i=1:Np
%                 for j=1:NTraces
%                     reconssig(i) = reconssig(i) + xa(i,j)*cos(vktheta(i,j))+xb(i,j)*sin(vktheta(i,j));
%                 end
%             end
            %     figure
            %     plot(xe(Nst:Nend)-reconssig(Nst:Nend));
            %     title('Differnce actual - real(Reconstructed)', 'FontSize', 11, 'Color', 'k')
            %     ylabel({'Diff.'}, 'FontSize', 10);

            figure(9)
            bsize=2;
            endplot=fix(Nend/bsize);
            bp=fix(2000/bsize);
            plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,1)),'-bs','LineWidth',1,'MarkerEdgeColor','b', 'MarkerSize',4,'MarkerIndices',1:bp:endplot);ylim([0 max(1.1*abs(FRFOrder(Nst:bsize:Nend,1,1)))]);
            hold on
            plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,2)),'-rd','LineWidth',1,'MarkerEdgeColor','r', 'MarkerSize',4,'MarkerIndices',1:bp:endplot)
            plot(blktime(Nst:bsize:Nend),abs(FRFOrder(Nst:bsize:Nend,3)),'-go','LineWidth',1,'MarkerEdgeColor','g', 'MarkerSize',4,'MarkerIndices',1:bp:endplot)
            plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,1)),':b','LineWidth',1.2);ylim([0 1.1*max((actamp(Nst:Nend,1)))]);
            plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,2)),':r','LineWidth',1.2);
            plot(blktime(Nst:Nend),abs(actamp(Nst:Nend,3)),':g','LineWidth',1.2);
            set(gcf,'position',[400,400,400,300])

            legend('1X',append(num2str(TraceOrder(2)),'X'),append(num2str(TraceOrder(3)),'Hz'),'actual');
            ylabel('Amp. (mm)');xlabel('time (s)')
            MAEa=sum(abs(acta(Nst:Nend,:)-xa(Nst:Nend,:)),1)/(Nend-Nst)
            MAEb=sum(abs(actb(Nst:Nend,:)-xb(Nst:Nend,:)),1)/(Nend-Nst)
            MPEA=100*max(abs(actamp(Nst:Nend,:)-abs(FRFOrder(Nst:Nend,:))))./max(actamp(Nst:Nend,:))
            mpmat(kk,:,opt)=MPEA;

        end
    end
end
% deltaf=0.2048624/rm(1,1);
% fc=deltaf*fs/2
return
fc=(fs/2)*0.2048624./rmmat;
ftick=[0.25 0.5 1 2 5 10 20 40 65];
for j=1:3
    figure()
    plot((fc),mpmat(:,j,1),'-bs','LineWidth',1,'MarkerEdgeColor','b', 'MarkerSize',4);
    hold on
    plot((fc),mpmat(:,j,2),'-rd','LineWidth',1,'MarkerEdgeColor','r', 'MarkerSize',4);
    set(gcf,'position',[400,400,400,300])
    xlim([0.25 70]);
    set(gca, 'XScale', 'log')
    set(gca,'XTick',ftick)
    legend('VKOT','SVKOT');
    ylabel(append('MPAE',num2str(j)));xlabel('cutoff freq.')
end
figure()
plot((fc),sum(mpmat(:,:,1),2)/3,'-bs','LineWidth',1,'MarkerEdgeColor','b', 'MarkerSize',4);
hold on
plot((fc),sum(mpmat(:,:,2),2)/3,'-rd','LineWidth',1,'MarkerEdgeColor','r', 'MarkerSize',4);
set(gcf,'position',[400,400,400,300])
xlim([0.25 70]);
set(gca, 'XScale', 'log')
set(gca,'XTick',ftick)
legend('VKOT','SVKOT');
ylabel('Avg. MPAE');xlabel('cutoff freq.')

return
Nst=1;
Nend=Np;
 figure()
 bsize=200;
 bp=12;
 endplot=fix(Nend/bsize);
 plot(ts(Nst:bsize:Nend),(f1(Nst:bsize:Nend)),'-bs','LineWidth',1,'MarkerEdgeColor','b', 'MarkerSize',4,'MarkerIndices',1:bp:endplot);ylim([0 max(1.1*f2(Nst:bsize:Nend))]);
 hold on
 plot(ts(Nst:bsize:Nend),(f2(Nst:bsize:Nend)),'-rd','LineWidth',1,'MarkerEdgeColor','r', 'MarkerSize',4,'MarkerIndices',1:bp:endplot)
 plot(ts(Nst:bsize:Nend),(f3(Nst:bsize:Nend)),'-go','LineWidth',1,'MarkerEdgeColor','g', 'MarkerSize',4,'MarkerIndices',1:bp:endplot)

 set(gcf,'position',[400,400,400,300])

 legend('1X',append(num2str(TraceOrder(2)),'X'),append(num2str(TraceOrder(3)),'Hz'));
 ylabel('Freq. (Hz)');xlabel('time (s)')

 figure()
 plot(ts,xe);
 ylim([-3 3]);
 set(gcf,'position',[400,400,400,300]) 
 ylabel('Amp. (mm)');xlabel('time (s)')
