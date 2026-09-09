@echo off
setlocal enabledelayedexpansion

:: Ruta base del proyecto
set "nxsa_AI_home=%cd%"
set "dataset=%nxsa_AI_home%\dataset"
set "obsid_file=%nxsa_AI_home%\utils\obsids_list_test.txt"
set "tmp_dir=%nxsa_AI_home%\tmp"

if not exist "%dataset%" mkdir "%dataset%"
if not exist "%tmp_dir%" mkdir "%tmp_dir%"
cd /d "%dataset%" || exit /b 1

for /f "usebackq delims=" %%A in ("%obsid_file%") do (
    for /f "tokens=1" %%B in ("%%A") do set "obs=%%B"

    if not exist "%dataset%\!obs!" (
        mkdir "%dataset%\!obs!\pps" 2>nul
        cd /d "%dataset%\!obs!\pps" || continue

        echo *** Downloading data for obsid: !obs! ...

        :: Events FIT file for PN
        curl -L -O -J "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=!obs!&name=PIEVLI&level=PPS&instname=PN&extension=FTZ" > "%tmp_dir%\!obs!_events.log" 2>&1

        :: Regions for all cameras
        curl -L -O -J "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=!obs!&name=REGION&level=PPS&extension=ASC" > "%tmp_dir%\!obs!_regions.log" 2>&1

        :: Image FIT file for PN
        curl -L -O -J "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=!obs!&name=IMAGE_8&level=PPS&instname=PN&extension=FTZ" > "%tmp_dir%\!obs!_image.log" 2>&1

    ) else (
        echo *** Folder !obs! exists ...
    )
)

endlocal