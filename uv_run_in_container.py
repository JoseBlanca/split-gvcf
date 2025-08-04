#!/usr/bin/env -S uv run python

from pathlib import Path
from subprocess import run
import shutil
import sys

CONTAINER_TEMPLATE_NAME = "split_genome"
EXTRA_MOUNTS = [(Path("/home/jose/tmp/vcfs_per_sample/"), Path("/vcfs"))]

CONTAINER_VENV = Path(".venv_container").absolute()
CONTAINER_VENV.mkdir(exist_ok=True)
PROJECT_DIR_IN_CONTAINER = Path("/code")
EXCLUDED_PROJECT_DIRS = [".venv", ".pytest_cache", CONTAINER_VENV.name]

HOME_DIR = Path.home().absolute()
PROJECT_DIR = Path(__name__).parent.absolute()
ROOT_HOME_DIR_IN_CONTAINER = Path("/root/")
UV_SHARE = ".local/share/uv"
UV_CACHE = ".cache/uv"


python_command = sys.argv[1:]
if not python_command:
    python_command = ["pytest", "test"]

PROJECT_MOUNTS = [
    (path, PROJECT_DIR_IN_CONTAINER / path.name)
    for path in PROJECT_DIR.iterdir()
    if path.name not in EXCLUDED_PROJECT_DIRS
]
PROJECT_MOUNTS.append((CONTAINER_VENV, PROJECT_DIR_IN_CONTAINER / ".venv"))
PROJECT_MOUNTS.append((HOME_DIR / UV_SHARE, ROOT_HOME_DIR_IN_CONTAINER / UV_SHARE))
PROJECT_MOUNTS.append((HOME_DIR / UV_CACHE, ROOT_HOME_DIR_IN_CONTAINER / UV_CACHE))

UV_COMMAND = ["uv", "run"] + python_command


VENV_DIR = PROJECT_DIR / ".venv"
if VENV_DIR.exists():
    shutil.rmtree(VENV_DIR)

CMD = ["podman", "run", "-it", "--rm"]
CMD.extend(["-v", f"{PROJECT_DIR}:{PROJECT_DIR_IN_CONTAINER}"])
for MOUNTS in (PROJECT_MOUNTS, EXTRA_MOUNTS):
    for local_dir, container_dir in MOUNTS:
        CMD.extend(["-v", f"{local_dir}:{container_dir}"])
CMD.extend(["-w", str(PROJECT_DIR_IN_CONTAINER)])
CMD.append(CONTAINER_TEMPLATE_NAME)
CMD.extend(UV_COMMAND)

run(CMD, check=True)
