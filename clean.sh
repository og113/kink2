#!/bin/bash

find data -type f -empty -exec rm {} \;
rm temp/*
rm data/*_step_0*
