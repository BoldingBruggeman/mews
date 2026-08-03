@echo off
REM close existing image windows
rem powershell -NoProfile -command "Get-Process | Where-Object {$_.ProcessName -like 'Photos*'} | Stop-Process -Force"

rem Step 1: copy output.nc to output1.nc (overwrite without asking)
copy /Y output.nc output1.nc

rem Step 2: run the model
gotm_sp1_2_2.exe

"C:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" "C:\Users\shaja\OneDrive - IOLR\MEWS\R_codes\Model_6_plots.R" 


@echo 
pause