#!/usr/bin/env bash
cat $@ | sed -E -e "s/chrM/chrMT/g" -e "s/chr//g"
