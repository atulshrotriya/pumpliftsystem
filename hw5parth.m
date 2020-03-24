function hw5parth
M = 1500; D=0.15; Cd=0.6; L=9; g=9.81; u=0.02; rho=1000;
d=0.015; be=4*(10^8); Pm=4.25*(10^7); Qm=0.0177;
Ac=pi*0.25*(D^2);
Pc0=(M*g)/Ac;
Ap=pi*0.25*(d^2);
[t,x]=ode45 (@eqns, [0 1], [Pc0 0 0 0]);
Pc=x(:,1); Q=x(:,2); y=x(:,3); dy=x(:,4);
figure(1)
plot(t,Q,'r','linewidth',2)
xlabel('time, s','fontsize',18)
ylabel('flow rate Q, m^3/s','fontsize',18)
figure(2)
plot(t,y,'r','linewidth',2)
xlabel('time, s','fontsize',18)
ylabel('elevator displacement y, m','fontsize',18)
figure(3)
plot(t,dy,'r','linewidth',2)
xlabel('time, s','fontsize',18)
ylabel('elevator velocity y'', m/s','fontsize',18)
figure(4)
plot(t,Pc/(M*g/Ac),'r','linewidth',2)
xlabel('time, s','fontsize',18)
ylabel('cylinder pressure P_c, N/m^2','fontsize',18)

    function dx=eqns(t,x)
        dx=zeros(4,1);
        Pc=x(1); Q=x(2); y=x(3); dy=x(4);
        Vc=(pi*0.25)*(((D^2)*y)+((d^2)*L));
        Pp=(Qm-Q)*(Pm/Qm);
        Rn=(4*rho*Q)/(pi*u*d);
        Pv=Pp-((rho*abs(Q)*Q)/(2*(Cd^2)*(Ap^2)));
        dx(1)=(be*(Q-Ac*dy))/Vc;
        if Rn>1187.6
            Res=((0.2414*(rho^0.75)*(u^0.25)*L)/(d^4.75))*Q*abs(Q)^0.75;
        else
            Res=((128*u*L)/(pi*(d^4)))*Q;
        end
        dx(2)=(Ap*((Pv-Pc)+(rho*g*L)-Res))/(rho*L);
        dx(3)=x(4);
        dx(4)=((Pc*Ac)/M)-g;
        Pc
    end
end

