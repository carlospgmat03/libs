function[rhoPPT] = PPT(rho,ModesAB)

rhoPPT=ppt_matrix(ModesAB) * rho *  ppt_matrix(ModesAB);



end

