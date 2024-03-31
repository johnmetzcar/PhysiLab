
Selecting TLGL-leukemia model initial conditions

Death attractor:

{'A20': 1, 'Apoptosis': 1, 'BID': 1, 'BclxL': 0, 'CREB': 0, 'CTLA4': 'X', 'Caspase': 1, 'Ceramide': 1, 'Cytoskeleton_signaling': 1, 'DISC': 1, 'ERK': 1, 'FLIP': 0, 'FYN': 1, 'Fas': 1, 'FasL': 1, 'FasT': 1, 'GAP': 0, 'GPCR': 0, 'GRB2': 1, 'GZMB': 1, 'IAP': 0, 'IFNG': 0, 'IFNGT': 1, 'IL15': 1, 'IL2': 0, 'IL2RA': 0, 'IL2RAT': 0, 'IL2RB': 1, 'IL2RBT': 1, 'JAK': 1, 'LCK': 1, 'MCL1': 0, 'MEK': 1, 'NFAT': 1, 'NFKB': 1, 'P2': 1, 'P27': 1, 'PDGF': 0, 'PDGFR': 0, 'PI3K': 1, 'PLCG1': 1, 'Proliferation': 0, 'RANTES': 1, 'RAS': 1, 'S1P': 0, 'SMAD': 0, 'SOCS': 0, 'SPHK1': 0, 'STAT3': 1, 'Stimuli': 1, 'TBET': 1, 'TCR': 'X', 'TNF': 1, 'TPL2': 1, 'TRADD': 0, 'ZAP70': 0, 'sFas': 0}

Survival attractors:


{'A20': 1, 'Apoptosis': 0, 'BID': 0, 'BclxL': 0, 'CREB': 0, 'CTLA4': 'X', 'Caspase': 0, 'Ceramide': 0, 'Cytoskeleton_signaling': 1, 'DISC': 0, 'ERK': 1, 'FLIP': 1, 'FYN': 1, 'Fas': 0, 'FasL': 1, 'FasT': 1, 'GAP': 0, 'GPCR': 1, 'GRB2': 1, 'GZMB': 1, 'IAP': 1, 'IFNG': 0, 'IFNGT': 1, 'IL15': 1, 'IL2': 0, 'IL2RA': 0, 'IL2RAT': 0, 'IL2RB': 1, 'IL2RBT': 1, 'JAK': 1, 'LCK': 1, 'MCL1': 1, 'MEK': 1, 'NFAT': 1, 'NFKB': 1, 'P2': 1, 'P27': 1, 'PDGF': 0, 'PDGFR': 1, 'PI3K': 1, 'PLCG1': 1, 'Proliferation': 0, 'RANTES': 1, 'RAS': 1, 'S1P': 1, 'SMAD': 1, 'SOCS': 0, 'SPHK1': 1, 'STAT3': 1, 'Stimuli': 1, 'TBET': 1, 'TCR': 'X', 'TNF': 1, 'TPL2': 1, 'TRADD': 0, 'ZAP70': 0, 'sFas': 1}

{'A20': 1, 'Apoptosis': 0, 'BID': 0, 'BclxL': 0, 'CREB': 0, 'CTLA4': 'X', 'Caspase': 0, 'Ceramide': 0, 'Cytoskeleton_signaling': 1, 'DISC': 0, 'ERK': 1, 'FLIP': 1, 'FYN': 1, 'Fas': 0, 'FasL': 1, 'FasT': 1, 'GAP': 0, 'GPCR': 1, 'GRB2': 1, 'GZMB': 1, 'IAP': 1, 'IFNG': 0, 'IFNGT': 1, 'IL15': 1, 'IL2': 0, 'IL2RA': 0, 'IL2RAT': 0, 'IL2RB': 1, 'IL2RBT': 1, 'JAK': 1, 'LCK': 1, 'MCL1': 1, 'MEK': 1, 'NFAT': 1, 'NFKB': 1, 'P2': 0, 'P27': 1, 'PDGF': 0, 'PDGFR': 1, 'PI3K': 1, 'PLCG1': 1, 'Proliferation': 0, 'RANTES': 1, 'RAS': 1, 'S1P': 1, 'SMAD': 1, 'SOCS': 0, 'SPHK1': 1, 'STAT3': 1, 'Stimuli': 1, 'TBET': 1, 'TCR': 'X', 'TNF': 1, 'TPL2': 1, 'TRADD': 0, 'ZAP70': 0, 'sFas': 1}

Notes on survival attractors:
Fix all nodes except:
CTLA4 - free - they ocsillate within the attractor
TCR - free - oscillate within the attractor
P2 - free -- this differeniator the two attractors - we will do put them in the attractor randomaly (50/50). 

Work from:

https://github.com/jcrozum/pystablemotifs/blob/master/Examples%20and%20Tutorials/Attractor%20Candidate%20Evaluation%20Advanced%20Tutorial.ipynb noting that we are running the network in this context:

PDGF.istate = 0;
IL15.istate = 1;
Stimuli.istate = 1;
Stimuli2.istate = 0;
CD45.istate = 0;
TAX.istate = 0;