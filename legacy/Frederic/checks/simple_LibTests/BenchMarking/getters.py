#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Function to generate data for benchMarking using the 1WCM molecule
# Author: Frederic Bonnet
# Date: 05/08/2016
#----------------------------------------------------------------------------
#System tools
import fnmatch
import os
import re
import sys
import platform
import multiprocessing
from sys import exit
#Data analysis packages
import numpy as np
################################################################################
# Getters for the package                                                      #
################################################################################
