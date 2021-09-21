#!/bin/bash

sphinx-apidoc -o source/ ../bin/
make html
