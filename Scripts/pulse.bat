@echo off
setlocal enabledelayedexpansion

REM Define your variable (e.g., a list of filenames)
set "gas=co2 ch4 n2o"
set "rcp=RCP26 RCP45 RCP60 RCP85"
set "cool=0 5 10 20 40"
set "timepulse=5 10 20 30"
REM Loop over the variable and execute a command 
for %%a in (%gas%) do (
    for %%b in (%cool%) do (
        for %%c in (%rcp%) do (
            for %%d in (%timepulse%) do (
    echo Processing file: gas %%a, cooling rate %%b, rcp %%c
    gams FAIR.gms --initial_conditions=historical_run --experiment=srmpulse --gas=%%a --cooling_rate=%%b --rcp=%%c --pulse_time=%%d --srm=1
) ) ) )

endlocal