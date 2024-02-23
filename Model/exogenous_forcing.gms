set t_rcp /1765*2500/;
set thisttot(t_rcp,t);
thisttot(t_rcp,t) = yes$(t_rcp.val-1764 eq t.val);

set sources /'solar','volcanic_annual','tropoz','stratoz','cloud_tot','totaer_dir','ch4oxstrath2o','landuse'/;
set sources /'fgassum','mhalosum'/;

parameter forcing_rcp(*,*,*);

$gdxin "input/RCPs_consolidated.gdx"
$load forcing_rcp=Forcing
$gdxin

forcing_exogenous(t)= sum( (t_rcp,sources)$thisttot(t_rcp,t),forcing_rcp(t_rcp,'%rcp%',sources)) + 3*sum(t_rcp$thisttot(t_rcp,t),forcing_rcp(t_rcp,'%rcp%','bcsnow'));
