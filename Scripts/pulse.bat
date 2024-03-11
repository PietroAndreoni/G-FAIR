@echo off
setlocal enabledelayedexpansion

REM Define your variable (e.g., a list of filenames)
set "gas=co2 ch4 n2o"

REM Loop over the variable and execute a command (e.g., printing each filename)
for %%i in (%gas%) do (
    echo Processing file: %%i
    gams FAIR.gms --initial_conditions=historical_run --experiment_ghg=pulse --gas=%%i --sai=1
    gams FAIR.gms --initial_conditions=historical_run --experiment_ghg=pulse --gas=%%i 
)

endlocal