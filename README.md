# CCBR Snakemake Pipeline Cookiecutter
This is a dummy folder framework for CCBR snakemake workflows.
New workflows can be started using this repository as a template.

## Creating PAT for GH 
This is a prerequisite for the next step. You will need [gh cli](https://cli.github.com/) installed on your laptop or use `/data/CCBR_Pipeliner/db/PipeDB/bin/gh_1.7.0_linux_amd64/bin/gh` on biowulf. Skip if can access github in an automated way already.

Personal Access Token (PAT) is required to access GitHub (GH) without having to authenticate by other means (like password) every single time. You can create a PAT by going [here](https://github.com/settings/tokens). Then you can copy the PAT and save it into a file on biowulf (say `~/gh_token`). Next, you can run the following command to set everything up correctly on biowulf (or your laptop)
```
gh auth login --with-token < ~/git_token
```

## Creating new repository
You can use [gh cli](https://cli.github.com/) to
 * create a new repository under CCBR, and
 * copy over the template code from CCBR_SnakemakePipelineCookiecutter
with the following command
```
gh repo create CCBR/<reponame> \
--description "<repo description>" \
--public \
--template CCBR/CCBR_SnakemakePipelineCookiecutter \
--confirm
```
On biowulf, you may have to specify the full path of the `gh` executable is located here: `/data/CCBR_Pipeliner/db/PipeDB/bin/gh_1.7.0_linux_amd64/bin/gh`

Then you can clone a local copy of the new repository:
```
gh repo clone CCBR/<reponame>.git
```

If you drop the `CCBR/` from the `gh` command above, then the new repo is created under your username. The commands would then look like this:
```
gh repo create <reponame> \
--description "<repo description>" \
--public \
--template CCBR/CCBR_SnakemakePipelineCookiecutter \
--confirm

gh repo clone <your_github_handle>/<reponame>.git
```

You can change `--public` to `--private` in the above `gh` command to make the newly created repository private.
