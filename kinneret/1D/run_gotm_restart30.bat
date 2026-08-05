@echo off
rem Step 1: copy output.nc to output1.nc and 30 years spinup restart (overwrite without asking)
copy /Y restart30.nc restart.nc
copy /Y output.nc output1.nc

rem Step 2: run the model
gotm_sp1_2.exe

@echo 
@echo 
@echo 
pause