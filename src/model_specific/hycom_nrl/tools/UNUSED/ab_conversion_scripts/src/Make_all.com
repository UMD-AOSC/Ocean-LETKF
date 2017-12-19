#!/bin/csh
make read_archive
make write_archive
