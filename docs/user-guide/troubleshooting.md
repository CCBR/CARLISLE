# Troubleshooting
Recommended steps to troubleshoot the pipeline.

## 1.1 Email
Check your email for an email regarding pipeline failure. You will receive an email from slurm@biowulf.nih.gov with the subject: Slurm Job_id=[#] Name=CARLISLE Failed, Run time [time], FAILED, ExitCode 1

## 1.2 Restart the run
After addressing the issue, unlock the output directory, perform another dry-run and check the status of the pipeline, then resubmit to the cluster.
```
#unlock dir
bash ./data/CCBR_Pipeliner/Pipelines/CARLISLE/carlisle --runmode=unlock --workdir=/path/to/output/dir

#perform dry-run
bash ./data/CCBR_Pipeliner/Pipelines/CARLISLE/carlisle --runmode=dryrun --workdir=/path/to/output/dir

#submit to cluster
bash ./data/CCBR_Pipeliner/Pipelines/CARLISLE/carlisle --runmode=run --workdir=/path/to/output/dir
```