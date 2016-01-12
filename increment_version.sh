#!/bin/bash

echo -n "Enter new version: "
read vers
sed -i "s/VERSION=.*/VERSION=$vers/" makefile
sed -i "s/Version:.*/Version: $vers/" package/DESCRIPTION
