classdef plib 
    
    methods(Static)
        %plug in your functions:
        %function[OUTPUT] = functionname[INPUT]
            % statement;
        %end;
        
        %GRID CONSTRUCTOR SPACE
        %function that constructs the space grid:
        function[xx] = xgrid(x_0,x_end,hx)           %set: x_0 = 0; x_end = 500; hx = 5;
            xx = x_0:hx:x_end;
        end
        
        %GRID CONSTRUCTOR TIME
        %function that constructs the time grid:
        function[tt] = tgrid(t_0,t_end,ht)           % set: t_0 = 0; t_end = 150; ht = 0.5;
             tt = t_0:ht:t_end;
        end
        

        %coeff  CONSTRUCTORfunciton that constructs the coefficients needed in the main Matrix
        function [c_m1,c_0,c_p1,omega] = coeffconstr(hx,v,D,alpha)
            c_m1    = -v/(2*hx)+D/(hx^2); 
            c_p1    = v/(2*hx)+D/(hx^2);
            c_0     = 2*D/hx^2; 
            omega   = -2*D/hx^2 - 2*alpha*hx*(-v/(2*hx)+D/hx^2); 
        end
        
        %MATRIX CONSTRUCTOR
        %Constructs the main matrix A given the input size = L/hx + 1 
        %and the coefficents of the diagonal elements of A. Those are given by the "coeffconstr(hx,v,D,alpha) function"
        function [A] = matrixconstr(size,c_m1,c_0,c_p1,omega)
            %Creating the Matrix, which is needed for the ode23 / CN
            A = gallery('tridiag',size,c_m1,-c_0,c_p1);
            A(1,1) = omega; A(1,2) = c_0;       %Applying Neumann bc at x_0
            A(end,end-1) = 0; A(end,end) = 1;   %Apllying Dirichlet bc at x_N
    
        end
        
        
        %SOURCE TERM %returns psix: shape(psix)= (length(xx),1) before t_switch, psi: the source for all space and time
        function [psix, source] = source(xx,tt,psi_0,t_switch,a,delta_a)
            t_0 = tt(1);
            source= zeros(length(xx),length(tt));
            psix = zeros(length(xx),1);
            syms z 
            Q = piecewise(a-delta_a<=z<=a+delta_a,psi_0/(2*t_switch*delta_a)*(1+cos(pi*(z-a)/delta_a)),0) ;
            psix(1:end) = subs(Q,z,xx);
            
            syms b
            R = piecewise(t_0<=b<=t_switch,psix,zeros(length(xx),1));
            source(:,:) = subs(R,b,tt);
        end
        
        
        %DIF EQ SOLVER ODE23
        %Input parameters are the Matrix A, see matrixconstr. function
        %The full description of our A matrix
        function [udot] = diffeq(t,u,A,sourcevec,tswitch)
         %for all times smaller than t0 we have leakage
            if t<tswitch
                 psi = sourcevec;%costant leakage over time
         %for all t>t0 the leakage is 0.
            else
                  psi = zeros(length(u),1);
            end
            udot = A*u+psi;
        end
        
       %DIF EQ SOLVER CRANK-NICELSON
       %input parameters are Matrix A(as T), xgird size, tgrid size, 
       %full source matrix , time grid , time-step and IC solution
       function [U] = CN(A_inv,B,xsize,tsize,source,tt,ht,u0)
           %T = old "A"

           U = zeros(tsize,xsize);
           u_current = u0;
           for cnt =1:length(tt)
               if cnt == 1
                    udot = u_current;
               elseif cnt < length(tt)
                    udot = (A_inv*B)*u_current + (ht/2).*A_inv*(source(:,1+cnt)+ source(:,cnt));
               else
                    udot = (A_inv*B)*u_current + (ht/2).*A_inv*(source(:,cnt)+ source(:,cnt-1));
               end
               U(cnt,:) = udot;
               u_current = udot; 
           end
       end
       
       %Alternative option for the CN function: 
       %{
        function [U] = CN(A_inv,B,xsize,tsize,source,tt,h_t,u)
           source = [source,source(:,end)]; %To avoid index out of bounce
           U = zeros(xsize,tsize);
           U(:,1) = u(:,1); %set initial condition to solution matrix
           for cnt =2:length(tt)
               U(:,cnt) = (A_inv*B)*U(:,cnt-1) + (h_t/2).*A_inv*(source(:,cnt+1)+ source(:,cnt));
           end
           U = U.'; 
        end
        %}
       
       %Creating matrices for Crank - Nicolson's solution
       function [A_inv,B] = CNmatrix(T,xsize,ht)
           %T = A form above
           Id = speye(xsize);
           A = Id - (ht/2)*T;
           B = Id + (ht/2)*T;  
           A_inv = inv(A);
       end
       
       %Create suraface plots for the soltuions
       function SurfPlot(x,t,sol,txt)
           [X,T] = meshgrid(x,t);
           figure
           surf(X,T,sol,'edgecolor','none')
           text(250,75,1,txt);
           xlabel("x-Distance from rivermouth[km]","FontSize",12);
           ylabel("t-Time[h]","FontSize",12);
           zlabel("u(x,t)[Mol/km]")
           title("Concentration of pesticide along a river with a pollution source","FontSize",14);
           colormap(jet)
       end
       
   
        %ANIMATION
        %Void function which shows animation once it is called. 
        %Input is the space grid, the time grid, and the solution u; 
        %If one wants to record, then the input variable 'record == 1'. 
       function animation(xx,tt,u,record)
        figure
        f1 = animatedline; 
        f1.Color = 'red';
        F = [];
        ymin = 0; ymax = max(u(:)) + 1/2*mean(u(:));
        axis([0 500 ymin ymax]);
        grid on;
        hold on
        title('Pesticide concentration in space over time');
        ylabel('$u(x,t)$ in $\left(\frac{mol}{km}\right)$');
        xlabel('Space ($x$) in km');
        dim = [.15 .1 .0 .8];
        for j = 1:length(tt)
            addpoints(f1,xx,u(j,:));                    
            annotation('textbox',dim,'String',['Time (t): ', num2str(tt(j))],'FitBoxToText','on'); %for time tag
            drawnow
            F = [F,getframe(gcf)]; 
            clearpoints(f1); %clear old lines to only plot temperature progress
            delete(findall(gcf,'type','annotation'));          %for non overlapping time tag
        end  
        hold off
        if record == 1
            video = VideoWriter('Animation.avi');
            open(video); writeVideo(video,F); close(video)
        end
        end
     
     %COMPUTE U_INFTY and DEAD FISH_INFTY
     function[deadinfty,U_infty] = Uinfty(U,p0,rho,xx,hx)
        U_infty = max(U);
        number = 0;
        for k = 1:length(xx)        %number = elements in space with pesticide above threshold
          if U_infty(k) > p0        %(number * hx) is the affected space above the threshold p0
              number = number+1;    
          end
        end
        deadinfty = rho * hx * number;   %compute number of dead fish
    end

    %COMPUTE U_TOT and DEAD FISH_TOT
    function[deadtot,U_tot] = Utot(time,U,p1,rho,xx,hx)
        U_tot = trapz(time,U,1);
        number = 0;
        for k = 1:length(xx)       %number = elements in space with pesticide above threshold
          if U_tot(k) > p1         %(number * hx) is the affected space above the threshold p0
              number = number+1;    
          end
        end
        deadtot = rho * hx * number;   %compute number of dead fish
    end
         
    %PLOT U_INFTY with INPUT: xx,U_infty,p0,rho,deadinfty
    function U_inftyPlot(xx,U_infty,p0,rho,deadinfty)
    figure
    plot(xx,U_infty,xx,p0*ones(1,length(xx)));
    title(['Maximal pesticide concentration $U_{\infty}$; Number of dead fish = ',num2str(deadinfty),'; $\rho$ = ',num2str(rho)])
    xlabel('Space [km]');
    ylabel('Pesticide concentration [mol/km]');
    legend('$U_{\infty}$',['Threshold p0 = ',num2str(p0)]);
    end
    
    %PLOT U_TOT with INPUT: xx,U_tot,p1,rho,deadtot 
    function  U_totPlot(xx,U_tot,p1,rho,deadtot)
    figure
    plot(xx,U_tot,xx,p1*ones(1,length(xx)));
    title(['Total pesticide concentration over time $U_{tot}$; Number of dead fish = ',num2str(deadtot),'; $\rho$ = ',num2str(rho)' ])
    xlabel('Space [km]');
    ylabel('$U_{tot}$ [mol*h/km]');
    legend('$U_{tot}$',['Threshold p1 = ',num2str(p1)]);
    end 
         
     
     
     %COMPUTE U FOR DIFFERENT RIVER VELOCITIES (all space/ all time / different v)
     %INPUT: Start velocity, end velocity, velocity steps, time grid, size of time grid, time steps, xgrid, size of xgrid, spacesteps, D, alpha, source for all time (matrix), initial condition u0
    function[UV] = VSol(vstart,vend,hv,tt,tsize,ht,xx,xsize,hx,D,alpha,source,u0)
        vvec = vstart:hv:vend;
        UV = zeros(length(tt),length(xx),length(vvec));
        for k = 1:length(vvec)    
        [c_m1,c_0,c_p1,omega] = plib.coeffconstr(hx,vvec(k),D,alpha);
        T = plib.matrixconstr(xsize,c_m1,c_0,c_p1,omega);
        [A_inv,B] = plib.CNmatrix(T,xsize,ht);
        u = plib.CN(A_inv,B,xsize,tsize,source,tt,ht,u0);
        UV(:,:,k) = u(:,:); 
        end
    end
    
    %COMPUTE U_tot AND NUMBER OF DEAD FISH FOR DIFFERENT RIVER VELOCITIES
    %INPUT: start velocity, end velocity, velocitysteps, xgrid, spacesteps, timegrid, UV=FULL SOLUTION FOR DIFFERENT VELOCITIES ([UV]=VSol(...)), threshold p1, rho
    function[deadtotV,U_totV] = UtotV(vstart,vend,hv,xx,hx,tt,UV,p1,rho)
        vvec = vstart:hv:vend;
        deadtotV = zeros(1,length(vvec));
        U_totV = zeros(length(vvec),length(xx));
        for k = 1:length(vvec)
           [deadtot,U_tot] = plib.Utot(tt,UV(:,:,k),p1,rho,xx,hx);
           deadtotV(k) = deadtot;
           U_totV(k,:) = U_tot(:);
        end
    end
    
    %COMPUTE U_INFTY AND NUMBER OF DEAD FISH FOR DIFFERENT RIVER VELOCITIES
    %INPUT: start velocity, end velocity, velocitysteps, xgrid, spacesteps, UV= FUll SOL for different velocities, threshold, fish denisity;
    function[deadinftyV,U_inftyV] = UinftyV(vstart,vend,hv,xx,hx,UV,p0,rho)
    vvec = vstart:hv:vend;
    deadinftyV = zeros(1,length(vvec));
    U_inftyV = zeros(length(vvec),length(xx));
        for k = 1:length(vvec)
           [deadinfty,U_infty] = plib.Uinfty(UV(:,:,k),p0,rho,xx,hx);
           deadinftyV(k) = deadinfty;
           U_inftyV(k,:) = U_infty;
        end
    end    
     
     %ANIMATION: UTOT FOR DIFFERENT RIVER VELOCITIES   
     %PLOTDIM TOT: [0 500 0 60];
     %PLOTDIM INFTY: [[0 500 0 1.5]]
     function UVanimation(U_totV,deadtotV,vstart,vend,hv,xx,p1,rho,plotdim)
        vvec = vstart:hv:vend;
        figure
        f1 = animatedline;
        f1.Color = 'red';
        F = [];
        axis(plotdim);
        grid on;
        hold on
        plot(xx,p1*ones(1,length(xx)))
        title(['$U_{\infty}$: Maximal pesticide concentration over time; $\rho = $',num2str(rho)]);
        ylabel('$U_{\infty}$ [mol/km]');
        xlabel('Space [km]');
        legend('$U_{\infty}(v)$',['Threshold p0 = ',num2str(p1)]);
        for j = 1:length(vvec)
            addpoints(f1,xx,U_totV(j,1:end));                    
            annotation('textbox','String',['River velocity v = ', num2str(vvec(j)),'[km/h]; Dead fish: ',num2str(deadtotV(j))],'FitBoxToText','on'); %for time tag
            drawnow
            pause(.01);
            F = [F,getframe(gcf)]; 
            clearpoints(f1);                                   %clear old lines to only plot temperature progress
            delete(findall(gcf,'type','annotation'));          %for non overlapping time tag
        end  
        hold off
         %video = VideoWriter('Uinfty.avi');
         %open(video); writeVideo(video,F); close(video)
        %}
     end
    
    %PLOT TIME TO OCAN AS FUNCTION OF RIVER VELOCITY
    function OceanvecPlt(vstart,vend,hv,Otime)
        figure
        plot(vstart:hv:vend,Otime)
        grid on
        xlabel('River velocity [km/h]');
        ylabel('Time pesticide to reach river mouth [h]');
        title('Pesticide to ocean travel time');
    end


    %CALCULATE: TIME TO OCEAN FOR DIFFERENT VELOCITIES
    function [Otime] = Time2Ocean(vstart,vend,hv,UV,tt)
        eps = 0.00001;
        Otime = zeros(1,length(vstart:hv:vend));
        for k = 1:length(vstart:hv:vend)
            index = find(UV(:,1,k) >= 0.001+eps,1);
            if ~isempty(index) 
                Otime(k) = tt(index) ;  
            end
        end
    end     
    
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
        
    end
end