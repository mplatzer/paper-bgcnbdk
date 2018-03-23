
# Pareto/GGG MCMC Draws

The Pareto/GGG MCMC draws were stored in a public Amazon AWS S3 bucket. These are large files, and thus are not contained in GitHub. In order to sync these your disk, you first need to install AWS CLI from https://docs.aws.amazon.com/cli/latest/userguide/installing.html and then execute:

> aws s3 sync --exclude online_grocery/* s3://phd-mcmc pggg-draws
