from pathlib import Path
from subprocess import run
import shutil

python_command = ["pytest", "test/"]
UV_COMMAND = ["uv", "run"] + python_command

PROJECT_DIR = Path(__name__).parent.absolute()
VENV_DIR = PROJECT_DIR / ".venv"
if VENV_DIR.exists():
    shutil.rmtree(VENV_DIR)

try:
    PROJECT_DIR_IN_CONTAINER = Path("/code")
    CONTAINER_TEMPLATE_NAME = "split_genome"
    CMD = ["podman", "run", "-it", "--rm"]
    CMD.extend(["-v", f"{PROJECT_DIR}:{PROJECT_DIR_IN_CONTAINER}"])
    CMD.extend(["-v", "uv-cache:/root/.cache/uv"])
    CMD.extend(["-w", str(PROJECT_DIR_IN_CONTAINER)])
    CMD.append(CONTAINER_TEMPLATE_NAME)
    CMD.extend(UV_COMMAND)
    print(" ".join(CMD))
    run(CMD, check=True)
except:
    raise
finally:
    if VENV_DIR.exists():
        shutil.rmtree(VENV_DIR)
