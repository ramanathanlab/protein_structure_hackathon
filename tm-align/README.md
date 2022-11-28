# TM-align module

## Instructions on Polaris

First, let's update the default conda environment location to be located on the performant `/lus/eagle` filesytem:
Add these lines to your `~/.condarc` file, where `<project-id>` and `<username>` correspond to your project and account:
```
pkgs_dirs:
  - /lus/eagle/projects/<project-id>/<username>/conda/pkgs
envs_dirs:
  - /lus/eagle/projects/<project-id>/<username>/conda/envs
env_prompt: ({name})
```
The last line simplifies the conda path in your prompt.

Now install the software:
```bash
git clone https://github.com/ramanathanlab/protein_structure_hackathon.git
module load conda/2022-09-08
conda create -y -n tm-align python=3.9
conda activate tm-align
conda install -y -c bioconda tmalign
conda install -y -c conda-forge tqdm
```

To run the software:
