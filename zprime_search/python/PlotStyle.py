import ROOT

# histogram properties -----------------------------------------------------------------------------
class Hist:
    def __init__(self, _histName, _fileName, _legendName, _color):
        self.histName = _histName
        self.fileName = _fileName
        self.legendName = _legendName
        self.color = _color

# draw properties ----------------------------------------------------------------------------------
class Draw:
    def __init__(self, _name, _titleX, _minX, _maxX, _titleY, _minY, _maxY, _bins, _log):
        self.name = _name
        self.titleX = _titleX
        self.minX = _minX
        self.maxX = _maxX
        self.titleY = _titleY
        self.minY = _minY
        self.maxY = _maxY
        self.bins = _bins
        self.log = _log

# get draw properties ------------------------------------------------------------------------------
def GetDrawP(name):
    drawPs = [
              Draw('ele_ll_pt', 'p_{T}(ee)', 0.0, 900.0, '', 0.5, 9e6, 10, 'y'),
              Draw('ele_ll_rapidity', 'Y(ee)', -5.0, 5.0, '', 0.5, 9e6, 4, 'y'),
              Draw('ele_ll_mass', 'M(ee)', 40.0, 2000.0, '', 0.5, 9e6, 10, 'xy'),
              Draw('hf_ll_pt', 'p_{T}(ee)', 0.0, 200.0, '', 0.5, 9e4, 4, 'y'),
              Draw('hf_ll_rapidity', 'Y(ee)', -5.0, 5.0, '', 0.5, 9e4, 4, 'y'),
              Draw('hf_ll_mass', 'M(ee)', 40.0, 900.0, '', 0.5, 9e4, 10, 'xy'),
              Draw('hf_hfEle_et', 'E_{T}(HF)', 20.0, 200.0, '', 0.5, 9e4, 2, 'y'),
              Draw('hf_hfEle_eta', '#eta(HF)', -5.2, 5.2, '', 0.5, 9e4, 4, 'y'),
              Draw('hf_raw_hfEle_eLong3x3', 'E9(L)', 100.0, 1600.0, '', 0.5, 9e4, 40, 'y'),
              Draw('hf_raw_hfEle_eLong5x5', 'E25(L)', 100.0, 1600.0, '', 0.5, 9e4, 40, 'y'),
              Draw('hf_raw_hfEle_eL9eL25', 'Iso Cut', 0.93, 1.01, '', 0.5, 9e4, 1, 'y'),
              Draw('hf_raw_hfEle_eCore', 'E(C)', 0.0, 1400.0, '', 0.5, 9e4, 40, 'y'),
              Draw('hf_raw_hfEle_eShort3x3', 'E9(S)', 0.0, 800.0, '', 0.5, 9e4, 40, 'y'),
              Draw('hf_raw_hfEle_cut2D', 'Bkg Cut', 0.2, 1.01, '', 0.5, 9e4, 20, 'y')
    ]
    for dp in drawPs:
        if name == dp.name:
            return dp
    return Draw(name, '', -1e4, 1e4, '', 0.5, 1e7, 1, 'y')

# get canvas width ---------------------------------------------------------------------------------
def GetW():
    return 560

# get canvas height --------------------------------------------------------------------------------
def GetH(flag='N'):
    if flag is 'N':
        return 400
    elif flag is 'R':
        return 160
    elif flag is 'S':
        return 240

# canvas style -------------------------------------------------------------------------------------
def SetCanvas(canvas, flag='N'):
    ROOT.gROOT.SetStyle('Plain')
    if flag in ['N', 'S']:
        canvas.SetMargin(0.10, 0.03, 0.12, 0.04)
    elif flag is 'R':
        canvas.SetMargin(0.10, 0.03, 0.20, 0.04)
    canvas.SetTicks()

# axes style ---------------------------------------------------------------------------------------
def SetAxes(hist, flag='N'):
    if flag is 'N':
        hist.GetXaxis().SetTitleSize(0.05)
        hist.GetXaxis().SetTitleOffset(1.0)
        hist.GetXaxis().SetLabelSize(0.04)
        hist.GetXaxis().SetLabelOffset(0.003)
        hist.GetYaxis().SetTitleSize(0.05)
        hist.GetYaxis().SetTitleOffset(0.9)
        hist.GetYaxis().SetLabelSize(0.04)
        hist.GetYaxis().SetLabelOffset(0.003)
    elif flag is 'R':
        hist.GetXaxis().SetTitleSize(0.12)
        hist.GetXaxis().SetTitleOffset(0.7)
        hist.GetXaxis().SetLabelSize(0.12)
        hist.GetYaxis().SetTitleSize(0.12)
        hist.GetYaxis().SetTitleOffset(0.4)
        hist.GetYaxis().SetLabelSize(0.12)
        hist.GetYaxis().SetNdivisions(507)
    elif flag is 'S':
        hist.GetXaxis().SetTitleSize(0.08)
        hist.GetXaxis().SetTitleOffset(0.6)
        hist.GetXaxis().SetLabelSize(0.08)
        hist.GetYaxis().SetTitleSize(0.08)
        hist.GetYaxis().SetTitleOffset(0.4)
        hist.GetYaxis().SetLabelSize(0.08)
        hist.GetYaxis().SetNdivisions(507)

# legend style -------------------------------------------------------------------------------------
def SetLegend(legend, flag='N'):
    if flag is 'N':
        legend.SetTextSize(0.04)
    elif flag is 'S':
        legend.SetTextSize(0.08)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
