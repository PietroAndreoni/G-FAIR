set sources /'solar','volcanic_annual','stratoz','cloud_tot','totaer_dir','landuse'/;
set sources /'fgassum','mhalosum'/;

parameter forcing_rcp(*,*,*);

$gdxin "input/RCPs_consolidated.gdx"
$load forcing_rcp=Forcing
$gdxin

forcing_exogenous(t)= (sum( (t_rcp,sources)$thisttot(t_rcp,t),forcing_rcp(t_rcp,'%rcp%',sources)) + 3*sum(t_rcp$thisttot(t_rcp,t),forcing_rcp(t_rcp,'%rcp%','bcsnow'))) / tstep;

** Exogenous forcing from SRM is 0 by default
forcing_srm(t) = 0;