%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sets all NaN, Inf, and erroneously large values to the fill value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var_out=set_fill_values(var_out,missing_value)

var_out.native.value(isinf(var_out.native.value))=missing_value;
var_out.native.value(isnan(var_out.native.value))=missing_value;
var_out.native.value(abs(var_out.native.value)>missing_value)=...
               sign(var_out.native.value(abs(var_out.native.value)>missing_value))*missing_value;
var_out.info.missing_value=missing_value;
