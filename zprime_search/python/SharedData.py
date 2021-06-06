# integrated luminosity ----------------------------------------------------------------------------
def Lumi(tag):
    if tag == 'ele':
        return 19.788 # +/- 2.6%
    if tag == 'hf':
        return 19.735 # +/- 2.6%

# scale factor for MC samples ----------------------------------------------------------------------
def McSF(tag):
    if tag == 'ele':
        return 1.0
    if tag == 'hf':
        return 0.63 # +/- 0.26

# scale factor for QCD estimation ------------------------------------------------------------------
def QcdSF(tag):
    if tag == 'ele':
        return 2.07 # +/- 0.39
    if tag == 'hf':
        return 1.69 # +/- 0.13
