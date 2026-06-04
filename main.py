import json
import os
import subprocess
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

BASE_DIR = Path(__file__).resolve().parent
SOLVER_PATH = Path(os.environ.get("SOLVER_PATH", BASE_DIR / "solver"))

app = FastAPI(title="Mixed Boundary Test Problem Solver")


@app.get("/api/solve")
def solve(
    n: int = Query(100, ge=2, le=20000),
    theta: float = Query(0.5, ge=0.0, le=1.0),
    gamma: float = Query(0.5, ge=0.0, le=1.0),
):
    if not SOLVER_PATH.exists():
        raise HTTPException(status_code=500, detail="C++ solver executable was not found")

    try:
        completed = subprocess.run(
            [str(SOLVER_PATH), str(n), str(theta), str(gamma)],
            check=False,
            capture_output=True,
            text=True,
            timeout=20,
        )
    except subprocess.TimeoutExpired as exc:
        raise HTTPException(status_code=504, detail="C++ solver timed out") from exc

    raw_output = completed.stdout.strip() or completed.stderr.strip()
    try:
        payload = json.loads(raw_output)
    except json.JSONDecodeError as exc:
        raise HTTPException(status_code=500, detail=f"Invalid solver JSON: {raw_output}") from exc

    if completed.returncode != 0 or "error" in payload:
        raise HTTPException(status_code=400, detail=payload.get("error", "Solver failed"))
    return payload


@app.get("/")
def index():
    return FileResponse(BASE_DIR / "index.html")


if (BASE_DIR / "static").exists():
    app.mount("/static", StaticFiles(directory=BASE_DIR / "static"), name="static")
