function[out] = generalized_state(ds, alphas, phis, string_specification, varargin );

if (nargin == 2)
	out = generalized_state(ds(1:4), ds(5:8), ds(9:end), alphas);
	return;
end
if (strcmp(string_specification , 'Morroco')) % {{{ 
	gamma_1 = Squeezed(ds(1)) + alphas(1)*eye(2);
	gamma_2 = Squeezed(ds(2)) + alphas(2)*eye(2);
	gamma_3 = Squeezed(ds(3)) + alphas(3)*eye(2);
	gamma_4 = Squeezed(ds(4)) + alphas(4)*eye(2);
	State = DirectSum({gamma_1 gamma_2 gamma_3 gamma_4});
	
	
	Gate = WhichModePhaseGate(phis(1), 1, 4); State = Gate*State*Gate';
	Gate = WhichModePhaseGate(phis(2), 3, 4); State = Gate*State*Gate';
	Gate = BeamSplitter([1 2], 4);         State = Gate*State*Gate';
	Gate = BeamSplitter([3 4], 4);         State = Gate*State*Gate';
	Gate = WhichModePhaseGate(phis(3), 3, 4); State = Gate*State*Gate';
	Gate = BeamSplitter([2 3], 4);         State = Gate*State*Gate';
	out = State;
	return;
end % }}}
if (strcmp(string_specification, 'With Interaction 1 and 4')) % {{{
	gamma_1 = Squeezed(ds(1)) + alphas(1)*eye(2);
	gamma_2 = Squeezed(ds(2)) + alphas(2)*eye(2);
	gamma_3 = Squeezed(ds(3)) + alphas(3)*eye(2);
	gamma_4 = Squeezed(ds(4)) + alphas(4)*eye(2);
	State = DirectSum({gamma_1 gamma_2 gamma_3 gamma_4});
	
	
	Gate = WhichModePhaseGate(phis(1), 1, 4); State = Gate*State*Gate';
	Gate = WhichModePhaseGate(phis(2), 3, 4); State = Gate*State*Gate';
	Gate = BeamSplitter([1 2], 4);         State = Gate*State*Gate';
	Gate = BeamSplitter([3 4], 4);         State = Gate*State*Gate';
	Gate = WhichModePhaseGate(phis(3), 3, 4); State = Gate*State*Gate';
	Gate = BeamSplitter([2 3], 4);         State = Gate*State*Gate';
	Gate = WhichModePhaseGate(phis(4), 1, 4); State = Gate*State*Gate';
	Gate = BeamSplitter([1 4], 4);         State = Gate*State*Gate';
	out = State;
	return;
end %}}}
if (strcmp(string_specification, 'James3070')) % {{{
	out  = InitialStateProduct(ds, alphas);
	out = ApplyPhase(out, phis(1), 2);
	out = ApplyBeamSplitter(out, [1 2], 0); % (1)
	out = ApplyBeamSplitter(out, [3 4], 0); % (2)
	out = ApplyPhase(out, phis(2), 3);      % phi_3
	out = ApplyPhase(out, phis(3), 4);      % phi_4
	% NumberForm[ArcCos[Sqrt[.3]], 8]
	% 0.99115659 
	% out = ApplyBeamSplitter(out, [2 4], 0); % Here depends on James email
	out = ApplyBeamSplitter(out, [2 4], 0, acos(sqrt(.3))); % Here depends on James email
	out = ApplyBeamSplitter(out, [1 3], 0); % Here depends on James email
	return;
end % }}}
if (strcmp(string_specification, 'JamesJan2010')) % {{{
	out  = InitialStateProduct(ds, alphas);
	out = ApplyPhase(out, phis(1), 1);                      % phi_1
	out = ApplyBeamSplitter(out, [1 2], 0);                 % BS (1)
	out = AttenuationChannel(out, acos(sqrt(.99)), 2);      % A1
	out = AttenuationChannel(out, acos(sqrt(.99)), 2);      % A2
	out = AttenuationChannel(out, acos(sqrt(.99)), 1);      % A3
	out = ApplyBeamSplitter(out, [3 4], 0);                 % BS(2)
	out = ApplyPhase(out, phis(2), 4);                      % phi_2
	out = ApplyBeamSplitter(out, [2 4], 0);                 % BS(3)
	out = AttenuationChannel(out, acos(sqrt(.99)), 4);      % A4
	out = AttenuationChannel(out, acos(sqrt(.99)), 2);      % A5
	out = AttenuationChannel(out, acos(sqrt(.99)), 2);      % A6
	out = ApplyPhase(out, phis(3), 3);                      % phi_3
	out = ApplyBeamSplitter(out, [1 3], 0, acos(sqrt(.3))); % BS(4)
	out = AttenuationChannel(out, acos(sqrt(.99)), 3);      % A7
	out = AttenuationChannel(out, acos(sqrt(.7)), 1);      % A8
	out = AttenuationChannel(out, acos(sqrt(.99)), 1);      % A9
	return;
end % }}}
if (strcmp(string_specification, 'JamesDec20092')) % {{{
	out  = InitialStateProduct(ds, alphas);

	%out = AttenuationChannel(out, acos(sqrt(.99)), 4);
	out = ApplyPhase(out, phis(1), 2);

	out = ApplyBeamSplitter(out, [1 2], 0); % (1)

	out = ApplyBeamSplitter(out, [3 4], 0); % (2)

	%out = AttenuationChannel(out, acos(sqrt(.99)), 4);
	out = ApplyPhase(out, phis(2), 3);      % phi_3

	%out = AttenuationChannel(out, acos(sqrt(.7)), 4);
	out = ApplyPhase(out, phis(3), 4);      % phi_4
	out = ApplyBeamSplitter(out, [2 4], 0, acos(sqrt(.3))); % 
	out = ApplyBeamSplitter(out, [1 3], 0); % Here depends on James email
	return;
end % }}}
if (strcmp(string_specification, 'JamesDec2009')) % {{{
	out  = InitialStateProduct(ds, alphas);

	out = AttenuationChannel(out, acos(sqrt(.99)), 4);
	out = ApplyPhase(out, phis(1), 2);

	out = ApplyBeamSplitter(out, [1 2], 0); % (1)

	out = ApplyBeamSplitter(out, [3 4], 0); % (2)

	out = AttenuationChannel(out, acos(sqrt(.99)), 4);
	out = ApplyPhase(out, phis(2), 3);      % phi_3


	out = AttenuationChannel(out, acos(sqrt(.7)), 4);
	
	out = ApplyPhase(out, phis(3), 4);      % phi_4
	out = ApplyBeamSplitter(out, [2 4], 0, acos(sqrt(.3))); % 
	out = ApplyBeamSplitter(out, [1 3], 0); % Here depends on James email
	return;
end % }}}
if (strcmp(string_specification, 'James5050')) % {{{
	out  = InitialStateProduct(ds, alphas);
	out = ApplyPhase(out, phis(1), 2);
	out = ApplyBeamSplitter(out, [1 2], 0); % (1)
	out = ApplyBeamSplitter(out, [3 4], 0); % (2)
	out = ApplyPhase(out, phis(2), 3);      % phi_3
	out = ApplyPhase(out, phis(3), 4);      % phi_4
	% NumberForm[ArcCos[Sqrt[.3]], 8]
	% 0.99115659 
	% out = ApplyBeamSplitter(out, [2 4], 0); % Here depends on James email
	out = ApplyBeamSplitter(out, [2 4], 0); % Here depends on James email
	out = ApplyBeamSplitter(out, [1 3], 0); % Here depends on James email
	return;
end % }}}
if (strcmp(string_specification, 'Phase lockers James3070')) % {{{
	AttenuationTheta = acos(sqrt(.97));

	out  = InitialStateProduct(ds, alphas);
	out = ApplyPhase(out, phis(1), 2);
	out = ApplyBeamSplitter(out, [1 2], 0); % (1)
	out = ApplyBeamSplitter(out, [3 4], 0); % (2)
	
	out = AttenuationChannel(out, AttenuationTheta, [2]);

	out = ApplyPhase(out, phis(2), 3);      % phi_3
	out = ApplyPhase(out, phis(3), 4);      % phi_4
	% NumberForm[ArcCos[Sqrt[.3]], 8]
	% 0.99115659 
	% out = ApplyBeamSplitter(out, [2 4], 0); % Here depends on James email
	out = ApplyBeamSplitter(out, [2 4], 0, acos(sqrt(.3))); % Here depends on James email
	out = ApplyBeamSplitter(out, [1 3], 0); % Here depends on James email
	out = AttenuationChannel(out, AttenuationTheta, [1]);
	out = AttenuationChannel(out, AttenuationTheta, [3]);
	return;
end % }}}
if (strcmp(string_specification, 'James3')) % {{{
	out  = InitialStateProduct(ds, alphas);
	out = ApplyPhase(out, phis(1), 2);
	out = ApplyBeamSplitter(out, [1 2], 0);
	out = ApplyBeamSplitter(out, [3 4], 0);
	out = ApplyPhase(out, phis(2), 3);
	out = ApplyPhase(out, phis(3), 4);
	out = ApplyBeamSplitter(out, [2 4], 0);
	out = ApplyBeamSplitter(out, [1 3], 0);
	return;
end % }}}
if (strcmp(string_specification, 'James no usar!!')) % {{{
	out  = InitialStateProduct(ds, alphas);
	out = ApplyPhase(out, phis(1), 1);
	out = ApplyBeamSplitter(out, [1 2], 0);
	out = ApplyBeamSplitter(out, [3 4], 0);
	out = ApplyPhase(out, phis(2), 3);
	out = ApplyPhase(out, phis(3), 4);
	out = ApplyBeamSplitter(out, [2 4], 0);
	out = ApplyBeamSplitter(out, [1 3], 0);
	return;
end % }}}
if (strcmp(string_specification, 'James2 on usar!!!')) % {{{
	out  = InitialStateProduct(ds, alphas);
	out = ApplyPhase(out, 2*pi-phis(1), 2);
	out = ApplyBeamSplitter(out, [1 2], 0);
	out = ApplyBeamSplitter(out, [3 4], 0);
	out = ApplyPhase(out, phis(2)-phis(1), 3);
	out = ApplyPhase(out, phis(3)-phis(1), 4);
	out = ApplyBeamSplitter(out, [2 4], 0);
	out = ApplyBeamSplitter(out, [1 3], 0);
	return;
end % }}}
