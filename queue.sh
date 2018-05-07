#!/bin/bash

squeue -t pd -o '%8u %8a %8Q' -S -p | uniq -c
