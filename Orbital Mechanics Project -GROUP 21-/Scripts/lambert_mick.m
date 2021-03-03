function [A,P,E,ERROR,VI,VF,TPAR,THETA]=lambert_mick(RI,RF,TOF,MU)
  
% 	    SUBROUTINE LAMBERT(RI,RF,TOF,MU,LONGWAY,A,P,E,ERROR,VI,VF,
%      C                   TPAR,THETA)
% ************************************************************************
% *
% *     This subroutine is a Lambert Algorithm which given two radius
% *     vectors and the time to get from one to the other, it finds the
% *     orbit connecting the two.  It solves the problem using a new
% *     algorithm developed by R. Battin.  It solves the Lambert problem
% *     for all possible types of orbits (circles, ellipses, parabolas,
% *     and hyperbolas).  The only singularity is for the case of a
% *     transfer angle of 360 degrees, which is a rather obscure case.
% *     It computes the velocity vectors corresponding to the given radius
% *     vectors except for the case when the transfer angle is 180 degrees
% *     in which case the orbit plane is ambiguous (an infinite number of
% *     transfer orbits exist).
%
% *     The algorithm computes the semi-major axis, the parameter (semi-
% *     latus rectum), the eccentricity, and the velocity vectors.
% *
% *     NOTE:  Open file 6 prior to calling LAMBERT.  If an error occurs
% *            or the 360 or 180 degree transfer case is encountered,
% *            LAMBERT writes to unit 6.
% *
% *     INPUTS TO THE SUBROUTINE
% *          RI(3)   = A three element array containing the initial
% *                    position vector (distance unit)
% *          RF(3)   = A three element array containing the final
% *                    position vector (distance unit)
% *          TOF     = The transfer time, time of flight (time unit)
% *          MU      = Gravitational parameter of primary
% *                    (distance unit)**3/(time unit)**2
% *          LONGWAY = Logical variable defining whether transfer is
% *                    greater or less than pi radians.
% *                       .TRUE.    Transfer is greater than pi radians
% *                       .FALSE.   Transfer is less than pi radians
% *
% *     OUTPUTS FROM THE SUBROUTINE
% *          A       = Semi-major axis of the transfer orbit
% *                    (distance unit)
% *          P       = Semi-latus rectum of the transfer orbit
% *                    (distance unit)
% *          E       = Eccentricity of the transfer orbit
% *          ERROR   = Error flag
% *                       .FALSE.   No error
% *                       .TRUE.    Error, routine failed to converge
% *          VI(3)   = A three element array containing the initial
% *                    velocity vector (distance unit/time unit)
% *          VT(3)   = A three element array containing the final
% *                    velocity vector (distance unit/time unit)
% *          TPAR    = Parabolic flight time between RI and RF (time unit)
% *          THETA   = The transfer angle (radians)
% *
% *          NOTE: The semi-major axis, positions, times, & gravitational
% *                parameter must be in compatible units.
% *
% *     MISSION PLANNING SUBROUTINES AND FUNCTIONS CALLED
% *         ABV, CROSSP, DOTP, PI, QCK(PI)
% *
% *     PROGRAMMER:    Chris D'Souza
% *
% *     DATE:          January 20, 1989
% *
% *     VERIFIED BY:   Darrel Monroe, 10/25/90
% *
% ************************************************************************

TOL = 1e-14;

PIE = pi;
TWOPI=2*pi;

% C
% C     ***  Compute the vector magnitudes and various
% C     ***  cross and dot products
% C

RIM2   = dot(RI,RI);
RIM    = sqrt(RIM2);            %modulo pos in
RFM2   = dot(RF,RF);
RFM    = sqrt(RFM2);            %modulo pos fin
CTH    = dot(RI,RF)/(RIM*RFM);  %coseno dell'angolo compreso
CR     = cross(RI,RF);          %normale al piano orbitale
STH    = norm(CR)/(RIM*RFM);    %modulo del versore momento della quantità di moto seno dell'angolo compreso

%*** choose angle for up angular momentum ***

if CR(3) < 0 %componente z del versore h nel ramo negativo Accade sia se vado all'indietro, sia se il ramo supera i 180°
	STH = -STH; %assumo il seno negativo (vettore diretto verso il basso)
end
THETA  = qck(atan2(STH,CTH));
B1     = sign(STH); if STH == 0; B1 = 1; end; %il 180 lo considero positivo
			
% C
% C     ***  Compute the chord and the semi-perimeter
% C
C= sqrt(RIM2 + RFM2 - 2*RIM*RFM*CTH); %modulo della distanza tra i due vett
S= (RIM + RFM + C)/2;
BETA   = 2*asin(sqrt((S-C)/S));
PMIN   = TWOPI*sqrt(S^3/(8*MU));
TMIN   = PMIN*(PIE-BETA+sin(BETA))/(TWOPI);
LAMBDA = B1*sqrt((S-C)/S);

% C
% C     ***  Compute L carefully for transfer angles less than 5 degrees
% C

if THETA*180/PIE <= 5  %trasformazione in gradi (angolo tra i due vettori inferiore a 5)
   W   = atan((RFM/RIM)^.25) - PIE/4;
   R1  = (sin(THETA/4))^2;
   S1  = (tan(2*W))^2;
   L   = (R1+S1)/(R1+S1+cos(THETA/2));
else
   L   = ((1-LAMBDA)/(1+LAMBDA))^2;
end

M= 8*MU*TOF^2/(S^3*(1+LAMBDA)^6);
TPAR   = (sqrt(2/MU)/3)*(S^1.5-B1*(S-C)^1.5);
L1     = (1 - L)/2;

% C
% C     ***  Initialize values of y, n, and x
% C

Y= 1;
N= 0;
ERROR  = 0;

if (TOF-TPAR) <= 1e-3
   X0  = 0;
else
   X0  = L;
end

X= -1.e8;

% C
% C     ***  Begin iteration
% C

while abs(X0-X) >= TOL && N <= 200
   N   = N+1;
   X   = X0;
   ETA = X/(sqrt(1+X) + 1)^2;

% C
% C  ***  Compute x by means of an algorithm devised by
% C  ***  Gauticci for evaluating continued fractions by the
% C  ***  'Top Down' method
% C

   DELTA = 1;
   U     = 1;
   SIGMA = 1;
   M1    = 0;
	 
   while abs(U) > TOL && M1 <= 200
		M1    = M1+1;
		GAMMA = (M1 + 3)^2/(4*(M1+3)^2 - 1);
		DELTA = 1/(1 + GAMMA*ETA*DELTA);
		U     = U*(DELTA - 1);
		SIGMA = SIGMA + U;
   end
	 
   C1 = 8*(sqrt(1+X)+1)/(3+1/(5 + ETA + (9*ETA/7)*SIGMA));
	 
   if N == 1
		DENOM = (1 + 2*X + L)*(3*C1 + X*C1 +4*X);
		H1 = (L+X)^2*(C1 + 1 + 3*X)/DENOM;
		H2 = M*(C1+X-L)/DENOM;
   else
	 
   	QR = sqrt(L1^2 + M/Y^2);
  	XPLL = QR - L1;
    LP2XP1 = 2*QR;
   	DENOM = LP2XP1*(3*C1 + X*C1+4*X);
   	H1 = ((XPLL^2)*(C1 + 1 + 3*X))/DENOM;
    H2 = M*(C1+X-L)/DENOM;
   end
	 
   B = 27*H2/(4*(1+H1)^3);
   U = -B/(2*(sqrt(B+1)+1));
	 
% C
% C  ***  Compute the continued fraction expansion K(u)
% C  ***  by means of the 'Top Down' method
% C

   DELTA = 1;
   U0 = 1;
   SIGMA = 1;
   N1 = 0;
   while N1 < 200 && abs(U0) >= TOL
			if N1 == 0
		   	GAMMA = 4/27;
   			DELTA = 1/(1-GAMMA*U*DELTA);
   			U0 = U0*(DELTA - 1);
   			SIGMA = SIGMA + U0;
			else
   
	 			for I8 = 1:2
					if I8 == 1
   					GAMMA = 2*(3*N1+1)*(6*N1-1)/(9*(4*N1 - 1)*(4*N1+1));
					else
   					GAMMA = 2*(3*N1+2)*(6*N1+1)/(9*(4*N1 + 1)*(4*N1+3));
					end

					DELTA = 1/(1-GAMMA*U*DELTA);
					U0 = U0*(DELTA-1);
					SIGMA = SIGMA + U0;
				end
			end
			
			N1 = N1 + 1;
   end
   
	 KU = (SIGMA/3)^2;
   Y = ((1+H1)/3)*(2+sqrt(B+1)/(1-2*U*KU));
   X0 = sqrt(((1-L)/2)^2+M/Y^2)-(1+L)/2;
	 
end

CONST = M*S*(1+LAMBDA)^2;
A = CONST/(8*X0*Y^2);

% C
% C     ***  Compute the velocity vectors
% C

if abs(THETA - PIE) >= 0.01
   R11 = (1 + LAMBDA)^2/(4*TOF*LAMBDA);
   S11 = Y*(1 + X0);
   T11 = (M*S*(1+LAMBDA)^2)/S11;
	 
   for K57=1:3
			VI(K57) = -R11*(S11*(RI(K57)-RF(K57))-T11*RI(K57)/RIM);
			VF(K57) = -R11*(S11*(RI(K57)-RF(K57))+T11*RF(K57)/RFM);
   end

%revised
elseif abs(THETA - PIE) < 0.01
    THETA=THETA+0.01;
    R11 = (1 + LAMBDA)^2/(4*TOF*LAMBDA);
    S11 = Y*(1 + X0);
    T11 = (M*S*(1+LAMBDA)^2)/S11;
	 
    for K57=1:3
			VI(K57) = -R11*(S11*(RI(K57)-RF(K57))-T11*RI(K57)/RIM);
			VF(K57) = -R11*(S11*(RI(K57)-RF(K57))+T11*RF(K57)/RFM);
   end
   %end revision

else
   APR=sqrt((S-RIM)*(1-(S-C)/(2*A)));
   TPRIME=TMIN-TOF;
   BPR=sign(TPRIME)*sqrt((S-RFM)*(1-S/(2*A)));
   COEFF=sqrt(2*MU/(RIM*C));
   VUCPR1=COEFF*APR;
   VUCMR1=COEFF*BPR;
   disp(' Since the transfer angle is almost 180 degrees')
   disp(' The velocity vectors cannot be uniquely resolved')
   disp(' They can be resolved as follows')
   disp([' vi =',num2str(VUCPR1),'*uc+r1   +',num2str(VUCMR1),'*uc-r1'])

end

if N1 == 200 | N == 200
   ERROR = 1;
   disp(' Lambert algorithm has not converged')
end

P = (2*RIM*RFM*Y^2*(1+X0)^2*sin(THETA/2)^2)/CONST;
E = sqrt(1 - P/A);
[A,P,E,ERROR,VI,VF,TPAR,THETA];
