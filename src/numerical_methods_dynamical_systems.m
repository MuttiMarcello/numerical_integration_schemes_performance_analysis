% Modeling and Simulation of Aerospace Systems
% Academic Year 2023/2024
% Assignment #1
% Mutti Marcello 220252 10698636

% Set the default font size for axes labels and tick labels
set(0,'DefaultAxesFontSize',20);

% Set the default font size for legend
set(0,'DefaultLegendFontSize',20);

% Set the default linewidth
set(0,'DefaultLineLineWidth',3);

%% EX1
clearvars; close all; clc;

% Function f
f=@(x) [x(2)^2-x(1)-2; -x(1)^2+x(2)+10];

% Analytical Jacobian Jf
Jf=@(x) [-1, 2*x(2); -2*x(1) 1];

% Approximated Jacobian (Forwards Differences)
Jf_f=@(x,e) [(f(x+e.*[1; 0])-f(x))./e, (f(x+e.*[0; 1])-f(x))./e];

% Reference solutions to the zero-finding problem
p2=[-1 0 4 1 6];
x2r=roots(p2);
x1r=x2r.^2-2;
z=[x1r.'; x2r.'];

% Newton's method initial conditions
x0=zeros(2,length(x1r));
x0(1,:)=ceil(x1r);
x0(2,:)=ceil(x2r);

% FD step sizes
D=10.^(0:-1:-6);

% Algorithm preallocation
n_it_an=zeros(1,length(x2r))';    % number of iterations of numerical solution
EA=zeros(size(n_it_an));          % error of analytical solution
EC=zeros(length(x2r),length(D));  % error of constrained FD solution
E=zeros(size(EC));                % error of FD solution
N=zeros(size(EC));                % number of iterations of FD solution

% Algorithm's stopping criteria
AbsTol=1e-9;    % Absolute Tolerance
RelTol=1e-12;   % Relative Tolerance
it_max=100;     % Limit of number of iterations 

% Analytical Solution
fprintf('Analytical Solutions\n')
for i=1:length(x2r)

    % Iteration startup
    it=1;
    abs_err=1;
    rel_err=1;
    xx_it=x0(:,i);
    
    while it<=it_max && abs_err>AbsTol && rel_err>RelTol

        % Newton's method iteration solution
        x_it=xx_it(:,it)-Jf(xx_it(:,it))\f(xx_it(:,it));

        % Absolute and relative error computation
        abs_err=norm(xx_it(:,it)-x_it);
        rel_err=abs_err/norm(x_it);

        it=it+1;
        xx_it=[xx_it x_it];

        % Missed convergence check
        if it==it_max
                fprintf('Iterations number limit reached')
        end
    end
    
    % Number of iterations of numerical solution
    n_it_an(i)=it-1;

    % Absolute Error of analytical solution
    EA(i)=norm(z(:,i)-xx_it(:,end));

    fprintf('Convergence to [%.4f%+.4fj %.4f%+.4fj] in %.0f iterations, with AE=%.4e\n',real(xx_it(1,end)),imag(xx_it(1,end)),real(xx_it(2,end)),imag(xx_it(2,end)),n_it_an(i),EA(i))
end

% FD constrained and unconstrained Solution
for j=1:length(D)
    for i=1:length(x2r)

        % Iteration startup
        it=1;
        xx_it=x0(:,i);
        
        while it<=n_it_an(i)

            % Newton's method iteration solution
            x_it=xx_it(:,it)-Jf_f(xx_it(:,it),D(j))\f(xx_it(:,it));

            it=it+1;
            xx_it=[xx_it x_it];
        end

        % Absolute error of constrained FD solution
        EC(i,j)=norm(z(:,i)-xx_it(:,end));
    end
    
    for i=1:length(x2r)

        % Iteration startup
        it=1;
        abs_err=1;
        rel_err=1;
        xx_it=x0(:,i);
        
        while it<=it_max && abs_err>AbsTol && rel_err>RelTol

            % Newton's method iteration solution
            x_it=xx_it(:,it)-Jf_f(xx_it(:,it),D(j))\f(xx_it(:,it));

            % Absolute and relative error computation
            abs_err=norm(xx_it(:,it)-x_it);
            rel_err=abs_err/norm(x_it);

            it=it+1;
            xx_it=[xx_it x_it];

            % Missed convergence check
            if it==it_max
                fprintf('Iterations number limit reached')
            end
        end

        % Number of iterations of numerical solution
        N(i,j)=it-1;
        
        % Absolute error of FD solution
        E(i,j)=norm(z(:,i)-xx_it(:,end)); 
    end
end

fprintf('\nFE Approximation, constrained iterations\n')
fprintf('Convergence with eps=%.0e and AE=%.4e\n',D(end),EC(1,end))
fprintf('Convergence with eps=%.0e and AE=%.4e\n',D(end),EC(2,end))
fprintf('Convergence with eps=%.0e and AE=%.4e\n',D(end),EC(3,end))
fprintf('Convergence with eps=%.0e and AE=%.4e\n',D(end),EC(4,end))

fprintf('\nFE Approximation\n')
fprintf('Convergence in %.0f iterations, with eps=%.0e and AE=%.4e\n',N(1,end),D(end),E(1,end))
fprintf('Convergence in %.0f iterations, with eps=%.0e and AE=%.4e\n',N(2,end),D(end),E(2,end))
fprintf('Convergence in %.0f iterations, with eps=%.0e and AE=%.4e\n',N(3,end),D(end),E(3,end))
fprintf('Convergence in %.0f iterations, with eps=%.0e and AE=%.4e\n',N(4,end),D(end),E(4,end))

figure
loglog(D,EC)
hold on
loglog(D,EA.*ones(size(EC)),'k-.','LineWidth',0.5)
grid minor
xlabel('\Delta')
ylabel('AE')
ylim([1e-16 1e-1])
legend('z_1','z_2','z_3','z_4','location','northwest')

figure
subplot(1,2,1)
loglog(D,E)
hold on
loglog(D,EA.*ones(size(EC)),'k-.','LineWidth',0.5)
grid minor
xlabel('\Delta')
ylabel('AE')
ylim([2e-16 6e-10])
legend('z_1','z_2','z_3','z_4','location','northwest')
subplot(1,2,2)
semilogx(D,N)
hold on
loglog(D,n_it_an.*ones(size(EC)),'k-.','LineWidth',0.5)
grid minor
xlabel('\Delta')
ylabel('n_{it}')
legend('z_1','z_2','z_3','z_4','location','northwest')

%% EX2
clearvars; close all; clc;

% ODE right hand side
f=@(t,x) x-2.*t.^2+2;

% Initial Condition
x0=1;

% Integration time interval
t_int=[0 2];

% Integrator step size
h=[0.5 0.2 0.05 0.01];

% Number of experiments for T_cpu computation
n_exp=100;

% Analytical solution
x_sol=@(t) 2.*t.^2+4.*t-exp(t)+2;
xx_sol=x_sol(linspace(t_int(1),t_int(2),100));  % Analytical solution array

% RK2 integration for different step sizes
[tt_1, xx_RK2_1] = RK2(f, x0, t_int, h(1));
[tt_2, xx_RK2_2] = RK2(f, x0, t_int, h(2));
[tt_3, xx_RK2_3] = RK2(f, x0, t_int, h(3));
[tt_4, xx_RK2_4] = RK2(f, x0, t_int, h(4));

% RK2 maximum absolute error computation
err_RK2=zeros(size(h));
err_RK2(1)=abs(x_sol(tt_1(end))-xx_RK2_1(end));
err_RK2(2)=abs(x_sol(tt_2(end))-xx_RK2_2(end));
err_RK2(3)=abs(x_sol(tt_3(end))-xx_RK2_3(end));
err_RK2(4)=abs(x_sol(tt_4(end))-xx_RK2_4(end));

% RK2 T_cpu realizations and average
T_RK2=zeros(length(h),n_exp);
for i=1:length(h)
    for j=1:n_exp
        T_RK2(i,j)=timeit(@() RK2(f, x0, t_int, h(i)),2);
    end
end
T_RK2=mean(T_RK2,2);

figure
subplot(2,2,1)
plot(linspace(t_int(1),t_int(2),100),xx_sol,'b',tt_1,xx_RK2_1,'-.x')
grid minor
title('h=0.5')
ylabel('x(t)')
legend('x(t)','RK2',Location='northwest')

subplot(2,2,2)
plot(linspace(t_int(1),t_int(2),100),xx_sol,'b',tt_2,xx_RK2_2,'-.x')
grid minor
title('h=0.2')
ylabel('x(t)')
legend('x(t)','RK2',Location='northwest')

subplot(2,2,3)
plot(tt_1,abs(x_sol(tt_1)-xx_RK2_1),'bo-')
grid minor
ylabel('AE')
xlabel('t')
legend('AE_{RK2}',Location='northwest')
fprintf('RK2 final AE for h=%.2f: %.4e\n',h(1),err_RK2(1))

subplot(2,2,4)
plot(tt_2,abs(x_sol(tt_2)-xx_RK2_2),'bo-')
grid minor
ylabel('AE')
xlabel('t')
legend('AE_{RK2}',Location='northwest')
fprintf('RK2 final AE for h=%.2f: %.4e\n',h(2),err_RK2(2))

figure
subplot(2,2,1)
plot(linspace(t_int(1),t_int(2),100),xx_sol,'b',tt_3,xx_RK2_3,'-.x')
grid minor
title('h=0.05')
ylabel('x(t)')
legend('x(t)','RK2',Location='northwest')

subplot(2,2,2)
plot(linspace(t_int(1),t_int(2),100),xx_sol,'b',tt_4,xx_RK2_4,'-.x')
grid minor
title('h=0.01')
ylabel('x(t)')
legend('x(t)','RK2',Location='northwest')

subplot(2,2,3)
plot(tt_3,abs(x_sol(tt_3)-xx_RK2_3),'bo-')
grid minor
ylabel('AE')
xlabel('t')
legend('AE_{RK2}',Location='northwest')
fprintf('RK2 final AE for h=%.2f: %.4e\n',h(3),err_RK2(3))

subplot(2,2,4)
plot(tt_4,abs(x_sol(tt_4)-xx_RK2_4),'bo-')
grid minor
ylabel('AE')
xlabel('t')
legend('AE_{RK2}',Location='northwest')
fprintf('RK2 final AE for h=%.2f: %.4e\n\n',h(4),err_RK2(4))

% RK4 integration for different step sizes
[~, xx_RK4_1] = RK4(f, x0, t_int, h(1));
[~, xx_RK4_2] = RK4(f, x0, t_int, h(2));
[~, xx_RK4_3] = RK4(f, x0, t_int, h(3));
[~, xx_RK4_4] = RK4(f, x0, t_int, h(4));

% RK4 maximum absolute error computation
err_RK4=zeros(size(h));
err_RK4(1)=abs(x_sol(tt_1(end))-xx_RK4_1(end));
err_RK4(2)=abs(x_sol(tt_2(end))-xx_RK4_2(end));
err_RK4(3)=abs(x_sol(tt_3(end))-xx_RK4_3(end));
err_RK4(4)=abs(x_sol(tt_4(end))-xx_RK4_4(end));

% RK4 T_cpu realizations and average
T_RK4=zeros(length(h),n_exp);
for i=1:length(h)
    for j=1:n_exp
        T_RK4(i,j)=timeit(@() RK4(f, x0, t_int, h(i)),2);
    end
end
T_RK4=mean(T_RK4,2);

figure
subplot(2,2,1)
plot(linspace(t_int(1),t_int(2),100),xx_sol,'b',tt_1,xx_RK4_1,'-.x')
grid minor
title('h=0.5')
ylabel('x(t)')
legend('x(t)','RK4',Location='northwest')

subplot(2,2,2)
plot(linspace(t_int(1),t_int(2),100),xx_sol,'b',tt_2,xx_RK4_2,'-.x')
grid minor
title('h=0.2')
ylabel('x(t)')
legend('x(t)','RK4',Location='northwest')

subplot(2,2,3)
plot(tt_1,abs(x_sol(tt_1)-xx_RK4_1),'bo-')
grid minor
ylabel('AE')
xlabel('t')
legend('AE_{RK4}',Location='northwest')
fprintf('RK4 final AE for h=%.2f: %.4e\n',h(1),err_RK4(1))

subplot(2,2,4)
plot(tt_2,abs(x_sol(tt_2)-xx_RK4_2),'bo-')
grid minor
ylabel('AE')
xlabel('t')
legend('AE_{RK4}',Location='northwest')
fprintf('RK4 final AE for h=%.2f: %.4e\n',h(2),err_RK4(2))

figure
subplot(2,2,1)
plot(linspace(t_int(1),t_int(2),100),xx_sol,'b',tt_3,xx_RK4_3,'-.x')
grid minor
title('h=0.05')
ylabel('x(t)')
legend('x(t)','RK4',Location='northwest')

subplot(2,2,2)
plot(linspace(t_int(1),t_int(2),100),xx_sol,'b',tt_4,xx_RK4_4,'-.x')
grid minor
title('h=0.01')
ylabel('x(t)')
legend('x(t)','RK4',Location='northwest')

subplot(2,2,3)
plot(tt_3,abs(x_sol(tt_3)-xx_RK4_3),'bo-')
grid minor
ylabel('AE')
xlabel('t')
legend('AE_{RK4}',Location='northwest')
fprintf('RK4 final AE for h=%.2f: %.4e\n',h(3),err_RK4(3))

subplot(2,2,4)
plot(tt_4,abs(x_sol(tt_4)-xx_RK4_4),'bo-')
grid minor
ylabel('AE')
xlabel('t')
legend('AE_{RK4}',Location='northwest')
fprintf('RK4 final AE for h=%.2f: %.4e\n',h(4),err_RK4(4))

figure
subplot(1,2,1)
loglog(h,T_RK2,'-o')
hold on
loglog(h,T_RK4,'-o')
grid minor
xlabel('h')
ylabel('T_{CPU} [s]')
legend('RK2','RK4')

subplot(1,2,2)
loglog(err_RK2,T_RK2,'-o')
hold on
loglog(err_RK4,T_RK4,'-o')
grid minor
ylabel('T_{CPU} [s]')
xlabel('AE_{max}')
legend('RK2','RK4')

%% EX3
clearvars; close all; clc;

% Reference linear ODE right hand side
A=@(a) [0 1; -1 2*cos(a)];

% RK2 linear operator
F_RK2=@(a,h) eye(2)+h*A(a)+1/2*(h*A(a))^2;

% RK4 linear operator
F_RK4=@(a,h) eye(2)+h*A(a)+1/2*(h*A(a))^2+1/6*(h*A(a))^3+1/24*(h*A(a))^4;

% Algorithm settings
m=200;      % number of discretized alpha values
n=10;       % h0 sweep vector size
Tol=1e-9;   % solution rejection tolerance

% Discretized alpha interval
alpha=linspace(pi,0,m);

% RK2 stability region problem initialization
hh_rk2=[];
aa_rk2=[];

% Solution for alpha=pi
h0=fzero(@(h) max(abs(eig(F_RK2(pi,h))))-1,5);
fprintf('RK2, alpha=pi, h=%.4f\n',h0)

% Iterating on alpha from pi to 0
for i=1:m

    % h0 sweep vector
    h0v=linspace(h0,0,n);

    for j=1:n-1

        % f_SR zero-finding solution
        h=fzero(@(h) max(abs(eig(F_RK2(alpha(i),h))))-1,h0v(j)+1);

        % First non-zero solution is collected
        if j==1 && h>Tol
            hh_rk2=[hh_rk2 h];
            aa_rk2=[aa_rk2 alpha(i)];
            h0=hh_rk2(end);

        % Different non-zero solutions are collected
        elseif abs(hh_rk2(end)-h)>Tol && h>Tol
            hh_rk2=[hh_rk2 h];
            aa_rk2=[aa_rk2 alpha(i)];
            h0=hh_rk2(end-1);
        end
    end
end

% Solution in the origin is added
hh_rk2=[hh_rk2 0];
aa_rk2=[aa_rk2 0];
hl_rk2=[hh_rk2.*cos(aa_rk2); hh_rk2.*sin(aa_rk2)];

% Completion by symmetry
hl_rk2=[hl_rk2(1,:) hl_rk2(1,end-1:-1:1);
        hl_rk2(2,:) -hl_rk2(2,end-1:-1:1)];

% RK4 stability region problem initialization
hh_rk4=[];
aa_rk4=[];

% Solution for alpha=pi
h0=fzero(@(h) max(abs(eig(F_RK4(pi,h))))-1,5);
fprintf('RK4, alpha=pi, h=%.4f\n',h0)

% Iterating on alpha from pi to 0
for i=1:m

    % h0 sweep vector
    h0v=linspace(h0,0,n);

    for j=1:n-1

        % f_SR zero-finding solution
        h=fzero(@(h) max(abs(eig(F_RK4(alpha(i),h))))-1,h0v(j)+1);

        % First non-zero solution is collected
        if j==1 && h>Tol
            hh_rk4=[hh_rk4 h];
            aa_rk4=[aa_rk4 alpha(i)];
            h0=hh_rk4(end);

        % Different non-zero solutions are collected
        elseif abs(hh_rk4(end)-h)>Tol && h>Tol
            hh_rk4=[hh_rk4 h];
            aa_rk4=[aa_rk4 alpha(i)];
            h0=hh_rk4(end-1);
        end
    end
end

% Solution in the origin is added
hh_rk4=[hh_rk4 0];
aa_rk4=[aa_rk4 0];
hl_rk4=[hh_rk4.*cos(aa_rk4); hh_rk4.*sin(aa_rk4)];

% TSP ordering
for i=1:length(hh_rk4)-1

    % Reference distance d
    d=norm(hl_rk4(:,i)-hl_rk4(:,i+1));

    for j=i+1:length(hh_rk4)-1

        % If a point is found at distance < d they are switched
        dn=norm(hl_rk4(:,i)-hl_rk4(:,j+1));
        if dn<d
            temp=hl_rk4(:,j+1);
            hl_rk4(:,j+1)=hl_rk4(:,i+1);
            hl_rk4(:,i+1)=temp;

            % Reference distance d update
            d=norm(hl_rk4(:,i)-hl_rk4(:,i+1));
        end
    end
end

% Completion by symmetry
hl_rk4=[hl_rk4(1,:) hl_rk4(1,end-1:-1:1);
        hl_rk4(2,:) -hl_rk4(2,end-1:-1:1)];

figure
plot(hl_rk2(1,:),hl_rk2(2,:))
hold on
grid minor
plot(hl_rk4(1,:),hl_rk4(2,:))
plot([-4 2],[0 0],'k-.',[0 0],[-3 3],'k-.','linewidth',0.5)
axis equal
xlim([-4 2])
ylim([-3 3])
xlabel('Re\{\lambda h\}')
ylabel('Im\{\lambda h\}')
legend('RK2','RK4','location','northwest')

%% EX4
clearvars; close all; clc;

% Reference linear ODE right hand side
A=@(a) [0 1; -1 2*cos(a)];

% RK1 linear operator
F_RK1=@(a,h) eye(2)+h*A(a);

% RK2 linear operator
F_RK2=@(a,h) eye(2)+h*A(a)+1/2*(h*A(a))^2;

% RK4 linear operator
F_RK4=@(a,h) eye(2)+h*A(a)+1/2*(h*A(a))^2+1/6*(h*A(a))^3+1/24*(h*A(a))^4;

% Initial Condition
x0=[1 1]';

% Integration time interval
t_int=[0 1];

% Number of discretized alpha values
m=250;

% Discretized alpha interval
alpha=linspace(0,pi,m);

% Target accuracy array
Tol=10.^(-(3:6));

% Algorithm preallocation
H_RK1=zeros(length(Tol),m);     % h solutions
H_RK2=H_RK1;
H_RK4=H_RK1;
LH_RK1=H_RK1;                   % h*lambda solutions
LH_RK2=LH_RK1;
LH_RK4=LH_RK1;
NF_RK1=zeros(1,length(Tol));    % number of function evaluations
NF_RK2=zeros(1,length(Tol));
NF_RK4=zeros(1,length(Tol));

% fzero initial conditions for RK1, RK2, RK4 respectively
h01=10.^(-3:-1:-6);
h02=0.045./3.^(0:3);
h04=0.6./2.^(0:3);

figure

% Iterating on Tol
for i=1:length(Tol)

    % Iterating on alpha from 0 to pi
    for j=1:m

        % f_AR zero-finding solution
        [h,~,ex_flag]=fzero(@(h) norm(expm(A(alpha(j)))*x0-((F_RK1(alpha(j),h))^((t_int(2)-t_int(1))/h))*x0,'inf')-Tol(i),h01(i));
        
        % Convergence check
        if ex_flag<=0
            error('fzero error')
        end

        % Function evaluations at alpha=pi
        if j==m
            NF_RK1(i)=ceil((t_int(2)-t_int(1))/h);
        end

        % Collection of solutions h
        H_RK1(i,j)=h;

        % Collection of lambda*h
        LH_RK1(i,j)=H_RK1(i,j)*exp(1i*alpha(j));
    end

    % Plot settings
    subplot(ceil(length(Tol)/2),ceil(length(Tol)/2),i)
    plot(real(LH_RK1(i,:)),imag(LH_RK1(i,:)),'b-')
    hold on
    plot(real(LH_RK1(i,:)),-imag(LH_RK1(i,:)),'b-')
    plot([1.1*min(real(LH_RK1(i,:))) 1.1*max(real(LH_RK1(i,:)))],[0 0],'k-.',[0 0],[-1.1*max(imag(LH_RK1(i,:))) 1.1*max(imag(LH_RK1(i,:)))],'k-.','linewidth',0.5)
    grid minor
    xlim([1.1*min(real(LH_RK1(i,:))) 1.1*max(real(LH_RK1(i,:)))])
    ylim([-1.1*max(imag(LH_RK1(i,:))) 1.1*max(imag(LH_RK1(i,:)))])
    xlabel('Re\{\lambda h\}')
    ylabel('Im\{\lambda h\}')
    str=sprintf('Tol=%.0e',Tol(i));
    legend(str,'Location','east')
end

figure

% Iterating on Tol
for i=1:length(Tol)

    % Iterating on alpha from 0 to pi
    for j=1:m

        % f_AR zero-finding solution
        [h,~,ex_flag]=fzero(@(h) norm(expm(A(alpha(j)))*x0-((F_RK2(alpha(j),h))^((t_int(2)-t_int(1))/h))*x0,'inf')-Tol(i),h02(i));
        
        % Convergence check
        if ex_flag<=0
            error('fzero error')
        end

        % Function evaluations at alpha=pi
        if j==m
            NF_RK2(i)=2*ceil((t_int(2)-t_int(1))/h);
        end

        % Colletion of solutions h
        H_RK2(i,j)=h;

        % Collection of solutions lambda*h
        LH_RK2(i,j)=H_RK2(i,j)*exp(1i*alpha(j));
    end

    % Plot settings
    subplot(ceil(length(Tol)/2),ceil(length(Tol)/2),i)
    plot(real(LH_RK2(i,:)),imag(LH_RK2(i,:)),'b-')
    hold on
    plot(real(LH_RK2(i,:)),-imag(LH_RK2(i,:)),'b-')
    plot([1.1*min(real(LH_RK2(i,:))) 1.1*max(real(LH_RK2(i,:)))],[0 0],'k-.',[0 0],[-1.1*max(imag(LH_RK2(i,:))) 1.1*max(imag(LH_RK2(i,:)))],'k-.','linewidth',0.5)
    grid minor
    xlim([1.1*min(real(LH_RK2(i,:))) 1.1*max(real(LH_RK2(i,:)))])
    ylim([-1.1*max(imag(LH_RK2(i,:))) 1.1*max(imag(LH_RK2(i,:)))])
    xlabel('Re\{\lambda h\}')
    ylabel('Im\{\lambda h\}')
    str=sprintf('Tol=%.0e',Tol(i));
    legend(str,'Location','northeast')
end

figure

% Iterating on Tol
for i=1:length(Tol)

    % Iterating on alpha from 0 to pi
    for j=1:m

        % f_AR zero-finding solution
        [h,~,ex_flag]=fzero(@(h) norm(expm(A(alpha(j)))*x0-((F_RK4(alpha(j),h))^((t_int(2)-t_int(1))/h))*x0,'inf')-Tol(i),h04(i));
        
        % Convergence check
        if ex_flag<=0
            error('fzero error')
        end

        % Function evaluations at alpha=pi
        if j==m
            NF_RK4(i)=4*ceil((t_int(2)-t_int(1))/h);
        end

        % Collection of solutions h
        H_RK4(i,j)=h;

        % Collection of lambda*h
        LH_RK4(i,j)=H_RK4(i,j)*exp(1i*alpha(j));
    end

    % Plot settings
    subplot(ceil(length(Tol)/2),ceil(length(Tol)/2),i)
    plot(real(LH_RK4(i,:)),imag(LH_RK4(i,:)),'b-')
    hold on
    plot(real(LH_RK4(i,:)),-imag(LH_RK4(i,:)),'b-')
    plot([1.1*min(real(LH_RK4(i,:))) 1.1*max(real(LH_RK4(i,:)))],[0 0],'k-.',[0 0],[-1.1*max(imag(LH_RK4(i,:))) 1.1*max(imag(LH_RK4(i,:)))],'k-.','linewidth',0.5)
    grid minor
    xlim([1.1*min(real(LH_RK4(i,:))) 1.1*max(real(LH_RK4(i,:)))])
    ylim([-1.1*max(imag(LH_RK4(i,:))) 1.1*max(imag(LH_RK4(i,:)))])
    xlabel('Re\{\lambda h\}')
    ylabel('Im\{\lambda h\}')
    str=sprintf('Tol=%.0e',Tol(i));
    legend(str,'Location','northeast')
end

figure
loglog(Tol,NF_RK4,'o-')
hold on
loglog(Tol,NF_RK2,'o-')
loglog(Tol,NF_RK1,'o-')
grid minor
xlabel('Tol')
ylabel('N_{F}')
legend('RK4','RK2','RK1')

%% EX5
clearvars; close all; clc;

% Reference linear ODE right hand side
A=@(a) [0 1; -1 2*cos(a)];

% BI2_theta linear operator
F_BI2=@(a,h,th) (eye(2)+(th-1)*h*A(a)+1/2*((th-1)*h*A(a))^2)\(eye(2)+h*th*A(a)+1/2*(h*th*A(a))^2);

% Algorithm settings
m=200;      % number of discretized alpha values
n=10;       % h0 sweep vector size
Tol=1e-9;   % solution rejection tolerance

% BI2_0.4 stability region problem initialization
% Theta
TH=0.4;

% Discretized alpha interval
alpha=linspace(0,pi,m);

hh=[];
aa=[];

% Initial guess for alpha=pi
h0=10;

% Iterating on alpha from 0 to pi
for i=1:m

    % h0 sweep vector
    h0v=linspace(h0,0,n);

    for j=1:n-1

        % f_SR zero-finding solution
        h=fzero(@(h) max(abs(eig(F_BI2(alpha(i),h,TH))))-1,h0v(j)+1);

        % First non-zero solution is collected
        if j==1 && h>Tol
            hh=[hh h];
            aa=[aa alpha(i)];
            h0=hh(end);

        % Different non-zero solutions are collected
        elseif abs(hh(end)-h)>Tol && h>Tol
            hh=[hh h];
            aa=[aa alpha(i)];
            h0=hh(end-1);
        end
    end
end

% Solution in the origin is added
hh=[hh 0];
aa=[aa 0];
hl=[hh.*cos(aa); hh.*sin(aa)];

% Completion by symmetry
hl=[hl(1,:) hl(1,end-1:-1:1);
    hl(2,:) -hl(2,end-1:-1:1)];

figure
plot(hl(1,:),hl(2,:))
hold on
plot([-2 12],[0 0],'k-.',[0 0],[-6 6],'k-.','linewidth',0.5)
grid minor
axis equal
xlim([-2 12])
ylim([-6 6])
xlabel('Re\{\lambda h\}')
ylabel('Im\{\lambda h\}')
title('BI2_{0.4}')

figure

% BI2_0.4 stability region problem initialization
% Theta values
th=[0.1 0.3 0.7 0.9];

% Iterating on theta
for k=1:length(th)

    % Choice of appropriate discretized alpha interval
    if th(k)<0.5
        alpha=linspace(0,pi,m);
    else
        alpha=linspace(pi,0,m);
    end

    hh=[];
    aa=[];

    % Initial guess for initial alpha
    h0=10;

    % Iterating on alpha
    for i=1:m

        % h0 sweep vector
        h0v=linspace(h0,0,n);

        for j=1:n-1

            % f_SR zero-finding solution
            h=fzero(@(h) max(abs(eig(F_BI2(alpha(i),h,th(k)))))-1,h0v(j)+1);

            % First non-zero solution is collected
            if j==1 && h>Tol
                hh=[hh h];
                aa=[aa alpha(i)];
                h0=hh(end);

            % Different non-zero solutions are collected
            elseif abs(hh(end)-h)>Tol && h>Tol
                hh=[hh h];
                aa=[aa alpha(i)];
                h0=hh(end-1);
            end
        end
    end

    % Solution in the origin is added
    hh=[hh 0];
    aa=[aa 0];
    hl=[hh.*cos(aa); hh.*sin(aa)];

    % Completion by symmetry
    hl=[hl(1,:) hl(1,end-1:-1:1);
        hl(2,:) -hl(2,end-1:-1:1)];
    
    plot(hl(1,:),hl(2,:))
    hold on
end
plot([-5.5 5.5],[0 0],'k-.',[0 0],[-5.5 5.5],'k-.','linewidth',0.5)
grid minor
axis equal
xlim([-5.5 5.5])
ylim([-5.5 5.5])
xlabel('Re\{\lambda h\}')
ylabel('Im\{\lambda h\}')
title('BI2_{\theta}')
legend('\theta=0.1','\theta=0.3','\theta=0.7','\theta=0.9','location','northwest')

%% EX6
clearvars; close all; clc;

% linead ODE matrix
B=[-180.5 +219.5;
   +179.5 -220.5];

% Initial condition
x0=[1 1]';

% Analytical solution
x_sol=@(t) exp(B*t)*x0;

% Integration time interval
t_int=[0 5];

% Integrator step size
h=0.1;

% Analytical solution computation
m=1000;

tt_sol=linspace(t_int(1),t_int(2),m);   % Solution over t_int
xx_sol=zeros(2,m);  
for i=1:m
    xx_sol(:,i)=expm(B*(tt_sol(i)-t_int(1)))*x0;
end

tt_solz=linspace(t_int(1),0.2,m);       % Zoomed solution over [0 0.2]
xx_solz=zeros(2,m);
for i=1:m
    xx_solz(:,i)=expm(B*(tt_solz(i)-t_int(1)))*x0;
end

% IVP matrix eigenvalues
lam=eig(B);

% Transients settling times
tau=log(100)./abs(lam);

figure
subplot(1,2,1)
plot(tt_sol,xx_sol(1,:),'b',tt_sol,xx_sol(2,:),'r')
grid minor
ylabel('x_i(t)')
xlabel('t')
legend('x_1(t)','x_2(t)')
title('Analytical solution')
subplot(1,2,2)
plot(tt_solz,xx_solz(1,:),'b',tt_solz,xx_solz(2,:),'r')
grid minor
ylabel('x_i(t)')
xlabel('t')
title('Fast dynamic zoom')

% RK4 integration
[tt,xx_rk4]=RK4(@(t,x) B*x,x0,t_int,h);

% IEX4 integration
[~,xx_iex4]=IEX4(@(t,x) B*x,x0,t_int,h);

% For the purpose of error analysis, the analytical solution is again
% computed with same time discretization as RK4,IEX4
m=length(tt);
tt_sole=linspace(t_int(1),t_int(2),m);
xx_sole=zeros(2,m);
for i=1:m
    xx_sole(:,i)=expm(B*(tt_sole(i)-t_int(1)))*x0;
end

figure
subplot(1,2,1)
plot(tt_sol,xx_sol,'b',tt,xx_rk4,'-.rx')
grid minor
ylabel('x_i(t)')
xlabel('t')
legend('x_i(t)','','RK4')
title('RK4')

subplot(1,2,2)
plot(tt_solz,xx_solz,'b',tt(1:3),xx_rk4(:,1:3),'-.rx')
grid minor
ylabel('x_i(t)')
xlabel('t')
title('RK4, fast dynamic zoom')

figure
subplot(1,2,1)
plot(tt_sol,xx_sol,'b',tt,xx_iex4,'rx-.')
grid minor
ylabel('x_i(t)')
xlabel('t')
legend('x_i(t)','','IEX4')
title('IEX4')

subplot(1,2,2)
plot(tt_solz,xx_solz,'b',tt(1:3),xx_iex4(:,1:3),'-.rx')
grid minor
ylabel('x_i(t)')
xlabel('t')
title('IEX4, fast dynamic zoom')

figure
subplot(1,2,1)
semilogy(tt,abs(xx_sole-xx_rk4),'-o')
grid minor
% ylim([1e-8 2e-4])
ylabel('AE')
xlabel('t')
legend('AE_{1,RK4}','AE_{2,RK4}','Location','northwest')

subplot(1,2,2)
semilogy(tt,abs(xx_sole-xx_iex4),'-o')
grid minor
ylim([1e-8 2e-4])
ylabel('AE')
xlabel('t')
legend('AE_{1,IEX4}','AE_{2,IEX4}','Location','northeast')

% reference linear ODE right hand side
A=@(a) [0 1; -1 2*cos(a)];

% RK4 linear operator
F_RK4=@(a,h) eye(2)+h*A(a)+1/2*(h*A(a))^2+1/6*(h*A(a))^3+1/24*(h*A(a))^4;

% IEX4 linear operator
F_IEX4=@(a,h) -1/6*((eye(2)-h*A(a))^(-1))+4*((eye(2)-1/2*h*A(a))^(-2))-27/2*((eye(2)-1/3*h*A(a))^(-3))+32/3*((eye(2)-1/4*h*A(a))^(-4));

% Algorithm settings
m=200;      % number of discretized alpha values
n=10;       % h0 sweep vector size
Tol=1e-9;   % solution rejection tolerance

% Discretized alpha interval
alpha=linspace(pi,0,m);

% RK4 stability region problem initialization
hh_rk4=[];
aa_rk4=[];

% Solution for alpha=pi
h0=fzero(@(h) max(abs(eig(F_RK4(alpha(1),h))))-1,5);

% Iterating on alpha from pi to 0
for i=1:m

    % h0 sweep vector
    h0v=linspace(h0,0,n);

    for j=1:n-1

        % f_SR zero-finding solution
        hit=fzero(@(h) max(abs(eig(F_RK4(alpha(i),h))))-1,h0v(j)+1);

        % First non-zero solution is collected
        if j==1 && hit>Tol
            hh_rk4=[hh_rk4 hit];
            aa_rk4=[aa_rk4 alpha(i)];
            h0=hh_rk4(end);

        % Different non-zero solutitons are collected
        elseif abs(hh_rk4(end)-hit)>Tol && hit>Tol
            hh_rk4=[hh_rk4 hit];
            aa_rk4=[aa_rk4 alpha(i)];
            h0=hh_rk4(end-1);
        end
    end
end

% Solution in the origin is added
hh_rk4=[hh_rk4 0];
aa_rk4=[aa_rk4 0];
hl_rk4=[hh_rk4.*cos(aa_rk4); hh_rk4.*sin(aa_rk4)];

% TSP ordering
for i=1:length(hh_rk4)-1

    % Reference distance d
    d=norm(hl_rk4(:,i)-hl_rk4(:,i+1));

    for j=i+1:length(hh_rk4)-1

        % If a point is found at distance < d they are switched
        dn=norm(hl_rk4(:,i)-hl_rk4(:,j+1));
        if dn<d
            temp=hl_rk4(:,j+1);
            hl_rk4(:,j+1)=hl_rk4(:,i+1);
            hl_rk4(:,i+1)=temp;

            % Reference distance d update
            d=norm(hl_rk4(:,i)-hl_rk4(:,i+1));
        end
    end
end

% Completion by symmetry
hl_rk4=[hl_rk4(1,:) hl_rk4(1,end-1:-1:1);
        hl_rk4(2,:) -hl_rk4(2,end-1:-1:1)];

% Discretized alpha interval
alpha=linspace(0,pi,m);

% IEX4 stability region problem initialization
hh_iex4=[];
aa_iex4=[];

% Solution for alpha=pi
h0=fzero(@(h) max(abs(eig(F_IEX4(alpha(1),h))))-1,15);

% Iterating on alpha from 0 to pi
for i=1:m

    % h0 sweep vector
    h0v=linspace(h0,0,n);

    for j=1:n-1

        % f_SR zero-finding solution
        hit=fzero(@(h) max(abs(eig(F_IEX4(alpha(i),h))))-1,h0v(j)+1);

        % Firsst non-zero solution is collected
        if j==1 && hit>Tol
            hh_iex4=[hh_iex4 hit];
            aa_iex4=[aa_iex4 alpha(i)];
            h0=hh_iex4(end);

        % Different non-zero solutions are collected
        elseif abs(hh_iex4(end)-hit)>Tol && hit>Tol
            hh_iex4=[hh_iex4 hit];
            aa_iex4=[aa_iex4 alpha(i)];
            h0=hh_iex4(end-1);
        end
    end
end

% Solution in the origin is added
hh_iex4=[hh_iex4 0];
aa_iex4=[aa_iex4 0];
hl_iex4=[hh_iex4.*cos(aa_iex4); hh_iex4.*sin(aa_iex4)];

% Completion by symmetry
hl_iex4=[hl_iex4(1,:) hl_iex4(1,end-1:-1:1);
        hl_iex4(2,:) -hl_iex4(2,end-1:-1:1)];

figure
plot(h*real(lam(2)),imag(lam(2)),'kx')
hold on
text(h*real(lam(2))+0.5,imag(lam(2))+0.5,'\lambda_2 h','FontSize',20)
plot(h*real(lam(1)),imag(lam(1)),'kx')
text(h*real(lam(1))+0.5,imag(lam(1))+0.5,'\lambda_1 h','FontSize',20)
plot(hl_rk4(1,:),hl_rk4(2,:),'r')
plot(hl_iex4(1,:),hl_iex4(2,:),'b')
plot([-41 14],[0 0],'k-.',[0 0],[-8 8],'k-.','linewidth',0.5)
grid minor
xlim([-41 14])
ylim([-8 8])
xlabel('Re\{\lambda h\}')
ylabel('Im\{\lambda h\}')

% Minimum h for RK4 convergence
h_rk4_st=hh_rk4(1)/abs(lam(2));
fprintf('Minimum step size required for RK4 convergence: %.8e\n',h_rk4_st)

%% EX7
clearvars; close all; clc;

% ODE right hand side
f=@(t,x) [-5/2*(1+8*sin(t))*x(1);
          (1-x(1))*x(2)+x(1)];

% Initial condition
x0=[1 1]';

% Integration time interval
t_int=[0 3];

% Integrator step size
h=0.1;

% AB3 integration
[tt,xx_ab3]=AB3(@(t,x) f(t,x),x0,t_int,h);

% AM3 integration
[~,xx_am3]=AM3(@(t,x) f(t,x),x0,t_int,h);

% ABM3 integration
[~,xx_abm3]=ABM3(@(t,x) f(t,x),x0,t_int,h);

% BDF3 integration
[~,xx_bdf3]=BDF3(@(t,x) f(t,x),x0,t_int,h);

% Reference ode113 solution for computation of absolute error
opt=odeset('RelTol',1e-12,'AbsTol',1e-12);
[~,xx_113]=ode113(@(t,x) f(t,x),tt,x0,opt);
xx_113=xx_113';

figure
subplot(1,2,1)
plot(tt,xx_113(1,:),'b')
hold on
plot(tt,xx_ab3(1,:),'-o')
grid minor
ylim([-0.2 1.1])
xlabel('t')
ylabel('x_1(t)')
legend('ode113','AB3','location','northeast')

subplot(1,2,2)
plot(tt,xx_113(2,:),'b')
hold on
plot(tt,xx_ab3(2,:),'-o')
grid minor
ylim([0 20])
xlabel('t')
ylabel('x_2(t)')
legend('ode113','AB3','location','northwest')

figure
subplot(1,2,1)
plot(tt,xx_113(1,:),'b')
hold on
plot(tt,xx_am3(1,:),'-*')
plot(tt,xx_abm3(1,:),'-x')
plot(tt,xx_bdf3(1,:),'-s')
grid minor
xlabel('t')
ylabel('x_1(t)')
ylim([-0.2 1.1])
legend('ode113','AM3','ABM3','BDF3','location','northeast')

subplot(1,2,2)
plot(tt,xx_113(2,:),'b')
hold on
plot(tt,xx_am3(2,:),'-*')
plot(tt,xx_abm3(2,:),'-x')
plot(tt,xx_bdf3(2,:),'-s')
grid minor
xlabel('t')
ylabel('x_2(t)')
ylim([0 20])
legend('ode113','AM3','ABM3','BDF3','Location','northwest')

figure
semilogy(tt,abs(xx_113(1,:)-xx_am3(1,:)),'-*')
hold on
semilogy(tt,abs(xx_113(1,:)-xx_abm3(1,:)),'-x')
semilogy(tt,abs(xx_113(1,:)-xx_bdf3(1,:)),'-s')
grid minor
xlabel('t')
ylabel('AE_1')
legend('AE_{AM3}','AE_{ABM3}','AE_{BDF3}','Location','southwest')

figure
semilogy(tt,abs(xx_113(2,:)-xx_am3(2,:)),'-*')
hold on
semilogy(tt,abs(xx_113(2,:)-xx_abm3(2,:)),'-x')
semilogy(tt,abs(xx_113(2,:)-xx_bdf3(2,:)),'-s')
grid minor
xlabel('t')
ylabel('AE_2')
legend('AE_{AM3}','AE_{ABM3}','AE_{BDF3}','Location','northwest')

%% FUNCTIONS

function [tt, xx] = RK2(f, x0, t_int, h)
%     Computes RK2 integration
%     Example: [tt, xx] = RK2(f, x0, t_int, h)
%     INPUTS:
%         f     [nx1] ODE right hand side f(t,x)
%         x0    [nx1] initial boundaty condition
%         t_int [1x2] itegration time interval
%         h     [1x1] step size
%     OUTPUTS:
%         tt    [1xN] discretized time interval (N dependent on h)
%         xx    [nxN] propagated solution
    
    % Input check
    if nargin<4
        fprintf('Not enough input arguments')
    end
    
    % Number of iterations N
    a=t_int(1);
    b=t_int(2);
    N=(b-a)/h;

    % Discretized time interval computation
    if N==round(N)
        tt=linspace(a,b,1+N);
    else
        tt=linspace(a,b,1+round(N));
        h=(b-a)/round(N);
        fprintf('Step has been rounded up to closest submultiple\n')
    end
    
    % Initial condition verticality check
    [m,n]=size(x0);
    if m<n
        x0=x0';
    end
    
    % Algorithm initialization
    xx=zeros(max(m,n),length(tt));
    xx(:,1)=x0;

    % Integration scheme
    for i=2:length(tt)
        xp=xx(:,i-1)+h*f(tt(i-1),xx(:,i-1));
        xx(:,i)=xx(:,i-1)+0.5*h*(f(tt(i-1),xx(:,i-1))+f(tt(i),xp));
    end
end

function [tt, xx] = RK4(f, x0, t_int, h)
%     Computes RK4 integration
%     Example: [tt, xx] = RK4(f, x0, t_int, h)
%     INPUTS:
%         f     [nx1] ODE right hand side f(t,x)
%         x0    [nx1] initial boundaty condition
%         t_int [1x2] itegration time interval
%         h     [1x1] step size
%     OUTPUTS:
%         tt    [1xN] discretized time interval (N dependent on h)
%         xx    [nxN] propagated solution
    
    % Input check
    if nargin<4
        fprintf('Not enough input arguments')
    end
    
    % Number of iterations N
    a=t_int(1);
    b=t_int(2);
    N=(b-a)/h;

    % Discretized time interval computation
    if N==round(N)
        tt=linspace(a,b,1+N);
    else
        tt=linspace(a,b,1+round(N));
        h=(b-a)/round(N);
        fprintf('Step has been rounded up to closest submultiple\n')
    end
    
    % Initial condition verticality check
    [m,n]=size(x0);
    if m<n
        x0=x0';
    end
    
    % Algorithm initialization
    xx=zeros(max(m,n),length(tt));
    xx(:,1)=x0;

    % Integration scheme
    for i=2:length(tt)
        k1=f(tt(i-1),xx(:,i-1));
        k2=f(tt(i-1)+0.5*h,xx(:,i-1)+0.5*h*k1);
        k3=f(tt(i-1)+0.5*h,xx(:,i-1)+0.5*h*k2);
        k4=f(tt(i),xx(:,i-1)+h*k3);
        xx(:,i)=xx(:,i-1)+(k1+2*k2+2*k3+k4)*h/6;
    end
    
end

function [tt, xx] = IEX4(f, x0, t_int, h)
%     Computes IEX4 integration
%     Example: [tt, xx] = IEX4(f, x0, t_int, h)
%     INPUTS:
%         f     [nx1] ODE right hand side f(t,x)
%         x0    [nx1] initial boundaty condition
%         t_int [1x2] itegration time interval
%         h     [1x1] step size
%     OUTPUTS:
%         tt    [1xN] discretized time interval (N dependent on h)
%         xx    [nxN] propagated solution
    
    % Input check
    if nargin<4
        fprintf('Not enough input arguments')
    end
    
    % Number of iterations N
    a=t_int(1);
    b=t_int(2);
    N=(b-a)/h;

    % Discretized time interval computation
    if N==round(N)
        tt=linspace(a,b,1+N);
    else
        tt=linspace(a,b,1+round(N));
        h=(b-a)/round(N);
        fprintf('Step has been rounded up to closest submultiple\n')
    end
    
    % Initial condition verticality check
    [m,n]=size(x0);
    if m<n
        x0=x0';
    end
    
    % Algorithm initialization
    xx=zeros(max(m,n),length(tt));
    xx(:,1)=x0;

    % Fsolve output is suppressed
    opt=optimoptions('fsolve','Display','off');

    % Integration scheme
    for i=2:length(tt)

        k1=fsolve(@(k) xx(:,i-1)+h*f(tt(i),k)-k,xx(:,i-1),opt);

        k2a=fsolve(@(k) xx(:,i-1)+1/2*h*f(tt(i-1)+1/2*h,k)-k,xx(:,i-1),opt);
        k2=fsolve(@(k) k2a+1/2*h*f(tt(i),k)-k,k2a,opt);

        k3a=fsolve(@(k) xx(:,i-1)+1/3*h*f(tt(i-1)+1/3*h,k)-k,xx(:,i-1),opt);
        k3b=fsolve(@(k) k3a+1/3*h*f(tt(i-1)+2/3*h,k)-k,k3a,opt);
        k3=fsolve(@(k) k3b+1/3*h*f(tt(i),k)-k,k3b,opt);

        k4a=fsolve(@(k) xx(:,i-1)+1/4*h*f(tt(i-1)+1/4*h,k)-k,xx(:,i-1),opt);
        k4b=fsolve(@(k) k4a+1/4*h*f(tt(i-1)+2/4*h,k)-k,k4a,opt);
        k4c=fsolve(@(k) k4b+1/4*h*f(tt(i-1)+4/4*h,k)-k,k4b,opt);
        k4=fsolve(@(k) k4c+1/4*h*f(tt(i),k)-k,k4c,opt);

        xx(:,i)=-1/6*k1+4*k2-27/2*k3+32/3*k4;
        
    end
    
end

function [tt, xx] = AB3(f, x0, t_int, h)
%     Computes AB3 integration
%     Example: [tt, xx] = AB3(f, x0, t_int, h)
%     INPUTS:
%         f     [nx1] ODE right hand side f(t,x)
%         x0    [nx1] initial boundaty condition
%         t_int [1x2] itegration time interval
%         h     [1x1] step size
%     OUTPUTS:
%         tt    [1xN] discretized time interval (N dependent on h)
%         xx    [nxN] propagated solution
   
    % Input check
    if nargin<4
        fprintf('Not enough input arguments')
    end
    
    % Number of iterations N
    a=t_int(1);
    b=t_int(2);
    N=(b-a)/h;

    % Discretized time interval computation
    if N==round(N)
        tt=linspace(a,b,1+N);
    else
        tt=linspace(a,b,1+round(N));
        h=(b-a)/round(N);
        fprintf('Step has been rounded up to closest submultiple\n')
    end
    
    % Initial condition verticality check
    [m,n]=size(x0);
    if m<n
        x0=x0';
    end
    
    % Algorithm initialization
    xx=zeros(max(m,n),length(tt));
    xx(:,1)=x0;

    % Integration scheme
    for i=2:length(tt)

        % RK3 startup
        if i<=3
            k1=f(tt(i-1),xx(:,i-1));
            k2=f(tt(i-1)+1/2*h,xx(:,i-1)+1/2*h*k1);
            k3=f(tt(i),xx(:,i-1)+h*(-k1+2*k2));
            xx(:,i)=xx(:,i-1)+h*(1/6*k1+2/3*k2+1/6*k3);

        % AB3 propagation
        else
            xx(:,i)=xx(:,i-1)+h/12*(23*f(tt(i-1),xx(:,i-1))-16*f(tt(i-2),xx(:,i-2))+5*f(tt(i-3),xx(:,i-3)));
        end
    end    
end

function [tt, xx] = AM3(f, x0, t_int, h)
%     Computes AM3 integration
%     Example: [tt, xx] = AM3(f, x0, t_int, h)
%     INPUTS:
%         f     [nx1] ODE right hand side f(t,x)
%         x0    [nx1] initial boundaty condition
%         t_int [1x2] itegration time interval
%         h     [1x1] step size
%     OUTPUTS:
%         tt    [1xN] discretized time interval (N dependent on h)
%         xx    [nxN] propagated solution
   
    % Input check
    if nargin<4
        fprintf('Not enough input arguments')
    end
    
    % Number of iterations N
    a=t_int(1);
    b=t_int(2);
    N=(b-a)/h;

    % Discretized time interval computation
    if N==round(N)
        tt=linspace(a,b,1+N);
    else
        tt=linspace(a,b,1+round(N));
        h=(b-a)/round(N);
        fprintf('Step has been rounded up to closest submultiple\n')
    end
    
    % Initial condition verticality check
    [m,n]=size(x0);
    if m<n
        x0=x0';
    end
    
    % Algorithm initialization
    xx=zeros(max(m,n),length(tt));
    xx(:,1)=x0;

    % Fsolve output is suppressed
    opt=optimoptions('fsolve','Display','off');
    
    % Integration scheme
    for i=2:length(tt)

        % RK3 startup
        if i<=2
            k1=f(tt(i-1),xx(:,i-1));
            k2=f(tt(i-1)+1/2*h,xx(:,i-1)+1/2*h*k1);
            k3=f(tt(i),xx(:,i-1)+h*(-k1+2*k2));
            xx(:,i)=xx(:,i-1)+h*(1/6*k1+2/3*k2+1/6*k3);

        % AM3 propagation
        else
            [xx(:,i),~,ex_flag]=fsolve(@(x) xx(:,i-1)+h/12*(5*f(tt(i),x)+8*f(tt(i-1),xx(:,i-1))-f(tt(i-2),xx(:,i-2)))-x,xx(:,i-1),opt);
            if ex_flag<=0
                error('AM3 fsolve error')
            end
        end
    end    
end

function [tt, xx] = ABM3(f, x0, t_int, h)
%     Computes ABM3 integration
%     Example: [tt, xx] = ABM3(f, x0, t_int, h)
%     INPUTS:
%         f     [nx1] ODE right hand side f(t,x)
%         x0    [nx1] initial boundaty condition
%         t_int [1x2] itegration time interval
%         h     [1x1] step size
%     OUTPUTS:
%         tt    [1xN] discretized time interval (N dependent on h)
%         xx    [nxN] propagated solution
   
    % Input check
    if nargin<4
        fprintf('Not enough input arguments')
    end
    
    % Number of iterations N
    a=t_int(1);
    b=t_int(2);
    N=(b-a)/h;

    % Discretized time interval computation
    if N==round(N)
        tt=linspace(a,b,1+N);
    else
        tt=linspace(a,b,1+round(N));
        h=(b-a)/round(N);
        fprintf('Step has been rounded up to closest submultiple\n')
    end
    
    % Initial condition verticality check
    [m,n]=size(x0);
    if m<n
        x0=x0';
    end
    
    % Algorithm initialization
    xx=zeros(max(m,n),length(tt));
    xx(:,1)=x0;

    % Integration scheme
    for i=2:length(tt)

        % RK3 startup
        if i<=3
            k1=f(tt(i-1),xx(:,i-1));
            k2=f(tt(i-1)+1/2*h,xx(:,i-1)+1/2*h*k1);
            k3=f(tt(i),xx(:,i-1)+h*(-k1+2*k2));
            xx(:,i)=xx(:,i-1)+h*(1/6*k1+2/3*k2+1/6*k3);

        % ABM3 propagation
        else
            xp=xx(:,i-1)+h/12*(23*f(tt(i-1),xx(:,i-1))-16*f(tt(i-2),xx(:,i-2))+5*f(tt(i-3),xx(:,i-3)));
            xx(:,i)=xx(:,i-1)+h/12*(5*f(tt(i),xp)+8*f(tt(i-1),xx(:,i-1))-f(tt(i-2),xx(:,i-2)));
        end
    end    
end

function [tt, xx] = BDF3(f, x0, t_int, h)
%     Computes BDF3 integration
%     Example: [tt, xx] = BDF3(f, x0, t_int, h)
%     INPUTS:
%         f     [nx1] ODE right hand side f(t,x)
%         x0    [nx1] initial boundaty condition
%         t_int [1x2] itegration time interval
%         h     [1x1] step size
%     OUTPUTS:
%         tt    [1xN] discretized time interval (N dependent on h)
%         xx    [nxN] propagated solution
   
    % Input check
    if nargin<4
        fprintf('Not enough input arguments')
    end
    
    % Number of iterations N
    a=t_int(1);
    b=t_int(2);
    N=(b-a)/h;

    % Discretized time interval computation
    if N==round(N)
        tt=linspace(a,b,1+N);
    else
        tt=linspace(a,b,1+round(N));
        h=(b-a)/round(N);
        fprintf('Step has been rounded up to closest submultiple\n')
    end
    
    % Initial condition verticality check
    [m,n]=size(x0);
    if m<n
        x0=x0';
    end
    
    % Algorithm initialization
    xx=zeros(max(m,n),length(tt));
    xx(:,1)=x0;

    % Fsolve output is suppressed
    opt=optimoptions('fsolve','Display','off');
    
    % Integration scheme
    for i=2:length(tt)

        % RK3 startup
        if i<=3
            k1=f(tt(i-1),xx(:,i-1));
            k2=f(tt(i-1)+1/2*h,xx(:,i-1)+1/2*h*k1);
            k3=f(tt(i),xx(:,i-1)+h*(-k1+2*k2));
            xx(:,i)=xx(:,i-1)+h*(1/6*k1+2/3*k2+1/6*k3);

        % BDF3 propagation
        else
            [xx(:,i),~,ex_flag]=fsolve(@(x) 18/11*xx(:,i-1)-9/11*xx(:,i-2)+2/11*xx(:,i-3)+6/11*h*f(tt(i),x)-x,xx(:,i-1),opt);
            if ex_flag<=0
                error('BDF3 fsolve error')
            end
        end
    end    
end