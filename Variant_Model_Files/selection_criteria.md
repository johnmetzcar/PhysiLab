we take 10 interventions from each method:
- pystablemotifs: take all 10
- IBMFA: from node states leading to apoptosis=1, select the 10 with apoptosis lowest when the control node is in the opposite state
- cubewalkers: i) only consider the strongest edge for each target ii) remove interventions within 2 hops of apoptosis  iii) take the top 10 remaining

Addiitonal criteria:
- Remove interventions present in both pystablemotifs and IBMFA. 
    - This is:
        - PDGFR-0
        - SPHK1-0
        - S1P-0
    - These are the single node interventions from Stable Motifs. 
    - Removing them from IBMFA list, but they will be written up as being predicted by both. 