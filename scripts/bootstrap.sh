#!/bin/bash
################################################################
#
# Bootstrap simple python repo / module / package
#
# author: Ruslan Ovsyannikov, 2023-2024
# 
# LICENSE: MIT
#
################################################################

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SCRIPT_NAME=$0
SCRIPT_ARGS=("$@")
MODULE_ROOT=$( cd ${SCRIPT_DIR} && cd .. && pwd )

REINSTALL=false

# Function to display help message
show_help() {
    echo "Usage: bootstrap.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -r,--reinstall    Delete existing .venv before creating a new one."
    echo "  -h,--help         Show this help message and exit."
}

# Check for arguments
for arg in "$@"; do
  case $arg in
    -r|--reinstall) REINSTALL=true; shift ;;
    -h|--help) show_help; exit 0 ;;
    *) echo "Unknown option: $arg"; show_help; exit 1 ;;
  esac
done

# sudo apt install python3 python3-venv

# check if venv module available
python3 -c "import ensurepip, venv" 2>/dev/null || { 
    echo "ERROR: The venv module or its dependencies are not available."
    echo "Please ensure that the required packages are installed." 
    exit -1
}

# Deactivate and remove existing virtual environment if --reinstall is specified
if [[ "$REINSTALL" = true && -d ${MODULE_ROOT}/.venv ]]; then
    if [[ "$VIRTUAL_ENV" == "${MODULE_ROOT}/.venv" ]]; then
        deactivate || echo "WARNING: Failed to deactivate virtual environment";
    fi
    rm -rf ${MODULE_ROOT}/.venv
    # Rerun the script with the same arguments
    echo "Reinstalling virtual environment..."
    exec "$SCRIPT_NAME" "${SCRIPT_ARGS[@]}"
fi

# Create virtual environment, if needed
if [[ ! -d ${MODULE_ROOT}/.venv ]]; then
    python3 -m venv ${MODULE_ROOT}/.venv || exit -1
fi

# Activate venv
if [[ ! -f ${MODULE_ROOT}/.venv/bin/activate ]]; then
    echo "Virtual environment does not exist or is damaged!" >&2
    exit -1
else
    . ${MODULE_ROOT}/.venv/bin/activate
fi

# Update pip - needed for older systems
pip install --upgrade pip

# Install dependencies/requirements
if [[ -f ${MODULE_ROOT}/requirements.txt ]]; then
    pip install -r ${MODULE_ROOT}/requirements.txt
fi

# Install development environment requirements
if [[ -f ${MODULE_ROOT}/dev-requirements.txt ]]; then
    pip install -r ${MODULE_ROOT}/dev-requirements.txt
fi

# Install the repo as a local package
if [[ -f ${MODULE_ROOT}/pyproject.toml ]]; then
    pip install -e ${MODULE_ROOT}
fi

# Install all local packages
if [[ -d ${MODULE_ROOT}/packages ]]; then
    for pckg in ${MODULE_ROOT}/packages/*/; do
        if [[ -d $pckg ]]; then
            pip install -e $pckg
        fi
    done
fi
