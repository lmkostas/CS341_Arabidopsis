#! /usr/bin/env python

import os
import subprocess

def set_env():
    source_env("prod_env.sh")

def source_env(filename):
    if os.path.isfile(filename):
        command = 'env -i bash -c "source {file} && env"'.format(file=filename)
        for line in subprocess.check_output(command, shell=True).split("\n"):
            if line:
                key, value = line.split("=")
                os.environ[key]= value
                if 'postgres' in value:
                    ip = value.split('@')[-1]
                    print "Using postgres {}".format(ip)

