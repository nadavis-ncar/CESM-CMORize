i%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert units 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=convert_units(data,ps,a,b,conversion)

if strcmp(conversion,'number_density_to_molar_mixing_ratio')
   pres=calculate_pressure(ps,a,b);

   rho=pres./(287*data.native_secondary.value);
   n_air=(6.02e23/0.028)*rho;

   data.native.value=(100^3)*data.native.value./n_air;
elseif strcmp(conversion,'molar_mixing_ratio_to_DU')
   data.native.value=data.native.value*(287*273.15/1e5);
end
