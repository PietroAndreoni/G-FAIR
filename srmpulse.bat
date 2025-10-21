setlocal enabledelayedexpansion

REM Define your variable (e.g., a list of filenames)
set "gas=co2 ch4 n2o"
set "rcp=RCP3PD RCP45 RCP6 RCP85"
set "cool=5 10 20 40"
set "timepulse=2 20 50"
set "stop=2200 2500 2900"
REM Loop over the variable and execute commands 
for %%a in (%gas%) do (
    for %%c in (%rcp%) do (
        for %%d in (%timepulse%) do (

            REM Print the main processing info (no cooling rate yet)
            echo Processing: gas=%%a, rcp=%%c, pulse_time=%%d

            REM Run the two GAMS calls that do not need cooling params
            gams FAIR.gms --initial_conditions=historical_run --experiment=pulse --gas=%%a --rcp=%%c --pulse_time=%%d
            gams FAIR.gms --initial_conditions=historical_run --experiment=srm --gas=%%a --rcp=%%c --pulse_time=%%d

            REM Now loop over cooling rates and stop times
            for %%b in (%cool%) do (
                for %%e in (%stop%) do (

                    REM compute f = %%e + 100 (batch arithmetic)
                    set /a f=%%e + 100

                    REM Use delayed expansion to reference f (use !f!)
                    echo Running with cooling_rate=%%b start=%%e end=!f! 

                    gams FAIR.gms --initial_conditions=historical_run --experiment=srm --gas=%%a --cooling_rate=%%b --start_rampdown=%%e --end_rampdown=!f! --rcp=%%c --pulse_time=%%d --srm_exogenous=1
                
                )
            )
        )
    )
)

pause

endlocal