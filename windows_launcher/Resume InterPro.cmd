@echo off
setlocal
pushd "%~dp0" || exit /b 1
where py >nul 2>&1
if errorlevel 1 (
    python launcher.py --resume-interpro
) else (
    py -3 launcher.py --resume-interpro
)
popd
pause
endlocal
