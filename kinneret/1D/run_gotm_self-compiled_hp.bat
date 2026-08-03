@echo off
set PATH=C:\Users\shaja\.conda\envs\fabmws2025\Library\bin;%PATH%
rem Step 1: copy output.nc to output1.nc and set restart=false (spinup run)
!copy /Y gotm_rf.yaml gotm.yaml
copy /Y output.nc output1.nc
rem Step 2: run the model for spinup
gotm_sp1_2.exe
rem set restart=true and run the model again
!copy /Y gotm_rt.yaml gotm.yaml
gotm_sp1_2.exe

@echo 
@echo 
@echo 
pause