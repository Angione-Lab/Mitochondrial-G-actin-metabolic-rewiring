function [modelRPMI,basisMediumRPMI]  = changeMedia_batch(model)
%changeMedia_batch
%RPMI 1640 media composition
%https://www.thermofisher.com/uk/en/home/technical-resources/media-formulation.114.html
%https://github.com/LukFil/SanityChecked/blob/cca19bda54b94db56f86e91643a15fdd2e882b61/RPMImediumSimulation.m
%https://github.com/SysBioChalmers/EnzymeConstrained_humanModels/blob/master/ComplementaryScripts/Simulation/setHamsMedium.m

% RPMI medium composition. 
medium_composition={'MAR09061'
'MAR09066'
'MAR09062'
'MAR09070'
'MAR09065'
'MAR09063' 
'MAR09071'
'MAR09067'
'MAR09038'
'MAR09039'
'MAR09040'
'MAR09041'
'MAR09042'
'MAR09043'
'MAR08386' 
'MAR09068'
'MAR09069'
'MAR09044'
'MAR09045'
'MAR09064'
'MAR09046'
'MAR09158'
'MAR09109'
'MAR09083'
'MAR09145'
'MAR09146'
'MAR09378'
'MAR09144'
'MAR09143'
'MAR09159'
'MAR09361'
'MAR09082'
'MAR09096'
'MAR09081'
'MAR09078'
'MAR09077'
'MAR09072'
'MAR01378'
'MAR09358'
'MAR09035'
'MAR09167'
'MAR09133'
'MAR09423'
'MAR09351'
'MAR09028'
};

% Medium concentrations
met_Conc_mM=[0.1
1.15
0.15
0.379
0.208
2 
0.136
0.133
0.0968
0.382
0.382
0.274
0.101
0.0909
0.153
0.174
0.286
0.168
0.0245
0.129
0.171
0.00863
0.00082
0.0214
0.000524
0.00227
0.082
0.00485
0.000532
0.00297
0.194
0.424
0
5.33
23.81
127.26
5.63
11.11
0
0
0
1
0
0.00326
0.0073 
];


current_inf = 1000;
set_inf =1000;

cellConc = 2.17*1e6;
t= 48;
cellWeight = 3.645e-12;


%% Definition of basic medium (defines uptake from the medium, not captured by the medium composition, all with same constraints)
mediumCompounds = {'MAR09058'; 'MAR09079'; 'MAR09047'; 'MAR09078'; 'MAR11420'; 'MAR09048'; 'MAR09072'; 'MAR09074' };

mediumCompounds_lb = -100;

%% Determine constraints to simulate medium: 
[modelRPMI,basisMediumRPMI] = setMediumConstraints(model, set_inf, current_inf, medium_composition, met_Conc_mM, cellConc, t, cellWeight, mediumCompounds, mediumCompounds_lb);
end