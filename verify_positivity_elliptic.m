function verify_positivity_elliptic(n, vol_Omega, p, c, rho, max_u_minus, supp_u_minus, Dm, norm_pstvpart, lambda1)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This program verifies the positivity of a exact solution u of -Laplacian u = f(u) on a bounded domain Omega,
	% where f is supposed to be a polynomial in this code.
	% Here, c is the coefficient matrix representing f. For example,
	% c = [4 , 1;
	%      5 , 2;
	%	   6 , 2.5;
	%	   -7, 3]
	% represents f(t) = 4t + 5|t|t + 6(|t|^1.5)t - 7t^3.
	% Therefore, the number of columns is always 2.
	% In fact, terms with negative coefficient can be ignored.
	% We suppose that the existence of u is proved nearby its approximation ^u
	% locally in ||\nabra(u - ^u)||_{L^2} <= rho, where rho is a positive numbear.
	% Before use this function, one should estimate with respect to ^u:
	% 	When f has superlinear terms with positive coefficients, e.g., f(t) = u^3,
	% 		- [Upper bound] max{^u-(x) : x in Omega}, where ^u- := max{0, -^u}.
	% 		- [Upper bound] Volume of supp(u-).
	% 	When f has a linear term, e.g., f(t) = 100(u-u^3),
	%		- [Upper bound] Volume of {x in Omega : ^u(x) <= m} given small m>0.
	% 		- [Lower bound] Lp norm of ^u+ over Dm for a subcritical exponent p, where ^u+ := max{0, ^u}.
	% Possible methods for computing the bounds can be found in the paper body.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% The inputs are the following:
	% 	n: [Integer] Dimension of domain Omega;
	%    currently restricted to n = 2 or 3.
	% 	vol_Omega: [Upper bound] Volume of entire domain Omega.
	% 	p: [Real] Given arbitrary subcritical exponent, e.g., p=2.
	% 	c: [Real] Coeficient matrix above.
	% 	rho: [Upper bound] Error between u and ^u in the sense of norm ||\nabra(u - ^u)||_{L^2}.
	%%%%%
	% When f does not have superlinear terms, set the followings arbitrary negative numbers, e.g., -999.
	% 	max_u_minus: [Upper bound] max{^u-(x) : x in Omega}, where ^u- := max{0, -^u}.
	% 	supp_u_minus: [Upper bound] Volume of supp(u-).
	%%%%%
	% When f does not have a linear term, set the followings arbitrary negative numbers, e.g., -999.
	% 	Dm: [Upper bound] Volume of {^u <= m} given small m>0.
	% 	norm_pstvpart: [Lower bound] Lp norm of ^u+ over Dm, where ^u+ := max{0, ^u}.
	%   lambda1: [Lower bound, Optional only for p=2] The minimal eigenvalue of -Laplacian over Omega.
	%            A sharp lower bound or the exact value make it easier to successfully verify positivity.
	%%%%%
	% Note: Replacing the above inputs with intervals will work in most cases.
	%       However, in order to avoid rare errors, please avoid to input intervals.	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% witten 29/06/2020 K.Tanaka
	% updated 07/09/2020 K.Tanaka and T.Asai
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	c(c(:,1)<0,:) = [];
	coef_lin = c(c(:,2)==1,1);
	if isempty(coef_lin) == 1
		coef_lin = 0;
	end
	c(c(:,2)==1,:)=[];

	left = intval(0);
	if isempty(c) == 0 % When f has superlinear terms with positive coefficients.
		for i = 1 : size(c,1)
			Cp2 = embedding(c(i,2)+1, vol_Omega, n);
			if Dm >= 0
				Cp1 = embedding(c(i,2)+1, Dm, n);
			else % Dm is not required so that it is input as negative Cp1 <= Cp2;
				Cp1 = Cp2;
			end
			tmp1 = intval(c(i,1))*(Cp1^2);
			tmp2 = intval(max_u_minus) * intval(supp_u_minus)^(1/intval(c(i,2)+1)) + Cp2*rho;
			left = left + tmp1*(tmp2^(c(i,2)-1));
		end
	end
	
	if coef_lin == 0 %When f does not have a linear term.
		fprintf('left <= %e\n', left.sup)
		fprintf('right = 1\n')
		if left.sup < 1
			disp('The positivity of solution u is confirmed.')
			return
		else
			error(['We cannot prove the inequality in Theorem 6.3.' newline 'Please recalculate constants with higher precision.'])
		end
	end

	if lambda1 < 0 || p ~= 2
		Cp_Omega = embedding(p, vol_Omega, n);
	else
		Cp_Omega = 1/sqrt(intval(lambda1));
	end
	fprintf('norm_pstvpart >= %e\n', norm_pstvpart)
	fprintf('Cp_Omega*rho <= %e\n', sup(Cp_Omega*rho))		
	if norm_pstvpart >= sup(Cp_Omega*rho)
		if n == 2
			lambda1_u_minus = intval('18.1684145355') / Dm;
		elseif n == 3
			lambda1_u_minus = intval('25.6463452794') / (Dm^(intval(2)/3));
		else
			error(['The current version is applicable for n = 2 or 3.' newline 'For the case in which n >= 4, refer to the paper body.'])
        end
        
        lambda1_u_minus

		fprintf('left <= %e\n', left.sup)
		fprintf('right >= %e\n', inf(1 - coef_lin/lambda1_u_minus))
		if left.sup < inf(1 - coef_lin/lambda1_u_minus)
			disp('The positivity of solution u is confirmed.')
			return
		else
			error(['We cannot prove the inequality in Theorem 4.2.' newline 'Please recalculate constants with higher precision.'])
		end
	else
		error(['We cannot prove norm_pstvpart >= Cp_Omega * rho.' newline 'Please recalculate rho and/or norm_pstvpart for lager m.'])
	end
end

function Cp_Omega = embedding(p, vol_Omega, n)
	if p == 2
		if n == 2
			Cp_Omega = sqrt(vol_Omega / intval('18.1684145355'));
			return
		elseif n == 3
			Cp_Omega = sqrt((vol_Omega^(intval(2)/3)) / intval('25.6463452794'));
			return
		else
			error(['The current version is applicable only for n = 2 or 3.' newline 'For the case in which n >= 4, refer to the paper body.'])
		end
	end
    Cp = talenti(p, n);
    p = intval(p);
    q = n*p/(p+n);
    vol_Omega = intval(vol_Omega);
    Cp_Omega = vol_Omega^((2-q)/2/q)*Cp;
end


function Cp = talenti(p, n)	
	n = intval(n);
    p = intval(p);
    q = n*p/(p+n);
    Pi = intval('pi');
    Cp = Pi^(-1/2)*n^(-1/q)*((q-1)/(n-q))^(1-1/q)*(gamma(1+n/2)*gamma(n)/gamma(n/q)/gamma(1+n-n/q))^(1/n);
end

