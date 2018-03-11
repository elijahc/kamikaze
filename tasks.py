from invoke import task
from subprocess import call

@task
def cov(ctx):

    call(['py.test','--cov=kamikaze','tests'])
