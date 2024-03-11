@echo off
setlocal enabledelayedexpansion

REM Define your variable (e.g., a list of filenames)
set "tremoval=20 40 60 80 100 120 140 160 180 200 220 240 260 280 300"

REM Loop over the variable and execute a command (e.g., printing each filename)
for %%i in (%tremoval%) do (
    echo Processing file: %%i
    gams FAIR.gms --initial_conditions=historical_run --experiment_ghg=pulse --gas=co2 --tremoval=%%i --sai=1
    gams FAIR.gms --initial_conditions=historical_run --experiment_ghg=pulse --gas=co2 --tremoval=%%i 
)

endlocal