#!/bin/bash

ENV_NAME="Lymphly"
ENV_PATH="./$ENV_NAME"

if ! command -v python3 &> /dev/null; then
    echo "Python was not found. Install Python first."
    exit 1
fi
if ! command -v pip &> /dev/null; then
    echo "pip was not found. Install pip first."
    exit 1
fi
if ! dpkg -s python3.10-venv &> /dev/null; then
    echo "Error: The package 'python3.10-venv' is not installed." >&2
    echo "Please install it by running:" >&2
    echo "  sudo apt update && sudo apt install python3.10-venv" >&2
    exit 1
fi

echo "Installing ipykernel..."
pip install --user ipykernel

echo "Creating an environment: $ENV_PATH"
python3 -m venv "$ENV_PATH"

source "$ENV_PATH/bin/activate"

echo "Updating pip inside the environment..."
pip install --upgrade pip

echo "Installing dependencies from requirements.txt ..."
pip install -r requirements.txt

echo "Adding Jupyter Kernel: $ENV_NAME"
python3 -m ipykernel install --user --name="$ENV_NAME" --display-name "$ENV_NAME"

echo "Done. The environment '$ENV_NAME' has been added to Jupyter kernels."
