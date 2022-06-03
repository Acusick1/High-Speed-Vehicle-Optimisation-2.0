function table = prandtl_meyer_table(g)
%PRANDTL_MEYER_TABLE Generate Prandtl-Meyer table

M = [1:0.0001:100, 100.01:0.01:500];

M1 = M(:);
g = g(:)';

vMp = ((g + 1)/(g - 1)).^0.5 ...
    * atan((((g - 1)/(g + 1)).*(M1.^2 - 1)).^0.5) ...
    - atan(((M1.^2) - 1).^0.5);

table = [M1, vMp];