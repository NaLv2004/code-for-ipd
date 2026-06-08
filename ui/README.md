# MIMO 2D Simulation UI

Start from this folder:

```powershell
python -m pip install -r requirements.txt
.\start_ui.ps1
```

Open:

```text
http://127.0.0.1:8000
```

The UI reads active `const` values from `platform_1/setting.h`. Starting a run writes the submitted values back to `setting.h`, builds `Release|x64` with MSBuild, copies the generated executable into `ui/runs/<run-id>/`, prepares an isolated `ui/runs/<run-id>/work/` directory, then runs that snapshot while streaming stdout/stderr to the browser.
