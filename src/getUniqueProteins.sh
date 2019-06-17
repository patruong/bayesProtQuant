#!/bin/bash

cat 500-PSSS3-equ\ decoy_Report.xls | cut -f 4 | sort | uniq -c > proteins.csv


