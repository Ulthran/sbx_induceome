FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_induceome_env

COPY envs/sbx_induceome_env.yml ./

# Install environment
RUN conda env create --file sbx_induceome_env.yml --name sbx_induceome

ENV PATH="/opt/conda/envs/sbx_induceome/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "sbx_induceome", "/bin/bash", "-c"]

# Run
CMD "bash"