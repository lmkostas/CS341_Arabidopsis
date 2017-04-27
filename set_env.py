#! /usr/bin/env python

import os
import subprocess

def set_env():
    filename = "prod_env.sh"
    if os.path.isfile(filename):
        command = 'env -i bash -c "source {file} && env"'.format(file=filename)
        for line in subprocess.check_output(command, shell=True).split("\n"):
            if line:
                key, value = line.split("=")
                os.environ[key]= value
