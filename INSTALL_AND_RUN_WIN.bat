@echo off
echo ===============================================
echo    GUIdedRNA - One-Click Installation
echo ===============================================
echo.
echo This will install and run GUIdedRNA.
echo Process may take 5-15 minutes.
echo.
pause

:: Create logs directory
if not exist logs mkdir logs

:: Check if Docker is installed
echo [1/4] Checking Docker...
docker --version >nul 2>&1
if %errorLevel% neq 0 (
    echo Docker not found! Please install Docker Desktop first:
    echo https://www.docker.com/products/docker-desktop/
    echo.
    echo After installing Docker Desktop:
    echo 1. Start Docker Desktop
    echo 2. Wait for it to fully load
    echo 3. Run this script again
    pause
    exit /b 1
)

:: Check if Docker is running
echo [2/4] Checking if Docker is running...
docker info >nul 2>&1
if %errorLevel% neq 0 (
    echo Docker is installed but not running!
    echo Please start Docker Desktop and try again.
    pause
    exit /b 1
)
echo ✓ Docker is ready!

:: Build the application
echo [3/4] Building GUIdedRNA (this may take 10-15 minutes)...
docker build -t guidedrna:latest .
if %errorLevel% neq 0 (
    echo ❌ Build failed! Check Docker Desktop is running.
    pause
    exit /b 1
)
echo ✓ Build completed!

:: Start the application
echo [4/4] Starting GUIdedRNA...

:: Stop any existing container
docker stop guidedrna-app >nul 2>&1
docker rm guidedrna-app >nul 2>&1

:: Find available port
set PORT=3838
:find_port
netstat -an | find "LISTENING" | find ":%PORT%" >nul
if %errorLevel% equ 0 (
    set /a PORT+=1
    if !PORT! leq 3850 goto :find_port
)

:: Create data and output directories if they don't exist
if not exist data mkdir data
if not exist output mkdir output

:: Build the docker run command with only existing drives
set DOCKER_CMD=docker run -d -p %PORT%:3838 -v "%cd%\data:/data" -v "%cd%\output:/output" -v "%USERPROFILE%:/host_home"

:: Check for existing drives and add them to the command
if exist "C:\" set DOCKER_CMD=%DOCKER_CMD% -v "C:\:/host_drives/C"
if exist "D:\" set DOCKER_CMD=%DOCKER_CMD% -v "D:\:/host_drives/D"
if exist "E:\" set DOCKER_CMD=%DOCKER_CMD% -v "E:\:/host_drives/E"
if exist "F:\" set DOCKER_CMD=%DOCKER_CMD% -v "F:\:/host_drives/F"
if exist "G:\" set DOCKER_CMD=%DOCKER_CMD% -v "G:\:/host_drives/G"
if exist "H:\" set DOCKER_CMD=%DOCKER_CMD% -v "H:\:/host_drives/H"

:: Complete the command
set DOCKER_CMD=%DOCKER_CMD% --name guidedrna-app guidedrna:latest

:: Execute the command
%DOCKER_CMD%

if %errorLevel% neq 0 (
    echo ❌ Failed to start container!
    echo Check Docker Desktop logs for more details.
    pause
    exit /b 1
)

echo.
echo ✅ GUIdedRNA is now running!
echo 🌐 Open your browser to: http://localhost:%PORT%
echo.
echo 📁 Available drives should now be visible in the file browser
echo 🏠 Your user folder is mounted as 'Host Home'
echo 💾 Available Windows drives are mounted as drive letters
echo.
echo To stop GUIdedRNA: go into docker desktop and close guidedrna-app
echo.
echo Waiting 10 seconds for the application to fully start...
timeout /t 10 /nobreak >nul
echo.
echo If the browser doesn't open automatically, manually go to:
echo http://localhost:%PORT%
echo.

:: Optional: Open browser automatically
start http://localhost:%PORT%

pause
