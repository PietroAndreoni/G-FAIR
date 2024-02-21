set t_rcp /1765*2500/;
set thisttot(t_rcp,t);
thisttot(t_rcp,t) = yes$(t_rcp.val-1764 eq t.val);

set sources /'solar_rf','volcanic_annual_rf','tropoz_rf','stratoz_rf','cloud_tot_rf','totaer_dir_rf','ch4oxstrath2o_rf','bcsnow_rf','landuse_rf'/;
set sources /'fgassum_rf','mhalosum_rf'/;

parameter forcing_rcp(*,*,*);

$gdxin "input/RCPs_consolidated.gdx"
$load forcing_rcp=Forcing
$gdxin

forcing_exogenous(t)= sum((t_rcp,sources),forcing_rcp(t_rcp,'%rcp%',sources)$thisttot(t_rcp,t));
