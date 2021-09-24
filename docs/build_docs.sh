#!/bin/bash

# This is just a quick script to generate the API docs then make the html.

sphinx-apidoc -o . ../bin/
make html
