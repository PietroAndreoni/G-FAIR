setlocal enabledelayedexpansion

REM Define your variable (e.g., a list of filenames)
set "gas=co2 ch4"
set "termination=10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000"

REM Loop over the variable and execute commands 
for %%a in (%gas%) do (

    REM Print the main processing info (no cooling rate yet)
    echo Processing: gas=%%a

    REM Run the two GAMS calls that do not need cooling params
    gams FAIR.gms --initial_conditions=historical_run --experiment=srm --gas=%%a --rcp=RCP45 --pulse_time=2

    for %%t in (%termination%) do (
        
        gams FAIR.gms --initial_conditions=historical_run --experiment=srm --gas=%%a --rcp=RCP45 --pulse_time=2 --termination_time=%%t

    )
)

pause

endlocal