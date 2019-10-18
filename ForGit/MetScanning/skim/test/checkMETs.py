
import ROOT
from DataFormats.FWLite import Events, Handle
from PhysicsTools.PythonAnalysis import *

events = Events(['skim.root'])
#events = Events(['root://eoscms.cern.ch//store/relval/CMSSW_7_4_4/RelValTTbar_13/GEN-SIM-RECO/PU50ns_MCRUN2_74_V8_38Tbis-v1/00000/AA6C4A0E-1E09-E511-839C-0025905A6084.root'])
#events = Events(['~/eos/cms/store/relval/CMSSW_7_4_4/RelValTTbar_13/GEN-SIM-RECO/PU50ns_MCRUN2_74_V8_38Tbis-v1/00000/AA6C4A0E-1E09-E511-839C-0025905A6084.root'])

edmCollections = [ 
  {'name':'pfMet', 'label':("pfMet"), 'edmType':"vector<reco::PFMET>"}, 
  {'name':'pfChMet', 'label':("pfChMet"), 'edmType':"vector<reco::PFMET>"}, 
  {'name':'pfCaloMet', 'label':("pfCaloMet",'','SKIM'), 'edmType':"vector<reco::PFMET>"}, 
  {'name':'HcalNoiseSummary', 'label':("hcalnoise"), 'edmType':"HcalNoiseSummary"}, 
]
events.toBegin()

def hbheNoiseIsoFilter(products):
  return \
    products["HcalNoiseSummary"].numIsolatedNoiseChannels()<10 and \
    products["HcalNoiseSummary"].isolatedNoiseSumE()<50. and \
    products["HcalNoiseSummary"].isolatedNoiseSumEt()<25.

handles={v['name']:Handle(v['edmType']) for v in edmCollections}
for i in range(10):
  events.to(i)
  products = {}
  for v in edmCollections:
    events.getByLabel(v['label'],handles[v['name']])
    products[v['name']] =handles[v['name']].product()

  if hbheNoiseIsoFilter(products):
    print "pfMet %f trkMet %f pfCaloMet %f"%(products['pfMet'][0].pt(),products['pfChMet'][0].pt(),products['pfCaloMet'][0].pt())
  else:
    print "Skip-> noise!"
