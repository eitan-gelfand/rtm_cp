* Encoding: UTF-8.

DATASET ACTIVATE DataSet2.
GLM Dprime.Asian.biasm.15 Dprime.Asian.biasm.21 Dprime.Asian.biasm.27 Dprime.Asian.biasm.33 
    Dprime.Asian.biasm.39 Dprime.Asian.biasm.45 Dprime.Asian.biasp.15 Dprime.Asian.biasp.21 
    Dprime.Asian.biasp.27 Dprime.Asian.biasp.33 Dprime.Asian.biasp.39 Dprime.Asian.biasp.45 
    Dprime.Caucasian.biasm.15 Dprime.Caucasian.biasm.21 Dprime.Caucasian.biasm.27 
    Dprime.Caucasian.biasm.33 Dprime.Caucasian.biasm.39 Dprime.Caucasian.biasm.45 
    Dprime.Caucasian.biasp.15 Dprime.Caucasian.biasp.21 Dprime.Caucasian.biasp.27 
    Dprime.Caucasian.biasp.33 Dprime.Caucasian.biasp.39 Dprime.Caucasian.biasp.45 BY Group
  /WSFACTOR=Race 2 Polynomial Regression 2 Polynomial Range 6 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(Race*Regression*Group) TYPE=BAR ERRORBAR=NO MEANREFERENCE=NO
  /EMMEANS=TABLES(OVERALL) 
  /EMMEANS=TABLES(Group) COMPARE ADJ(BONFERRONI)
  /EMMEANS=TABLES(Race) COMPARE ADJ(BONFERRONI)
  /EMMEANS=TABLES(Regression) COMPARE ADJ(BONFERRONI)
  /EMMEANS=TABLES(Range) COMPARE ADJ(BONFERRONI)
  /EMMEANS=TABLES(Group*Race) 
  /EMMEANS=TABLES(Group*Regression) 
  /EMMEANS=TABLES(Group*Range) 
  /EMMEANS=TABLES(Race*Regression) 
  /EMMEANS=TABLES(Race*Range) 
  /EMMEANS=TABLES(Regression*Range) 
  /EMMEANS=TABLES(Group*Race*Regression) 
  /EMMEANS=TABLES(Group*Race*Range) 
  /EMMEANS=TABLES(Group*Regression*Range) 
  /EMMEANS=TABLES(Race*Regression*Range) 
  /EMMEANS=TABLES(Group*Race*Regression*Range) 
  /PRINT=DESCRIPTIVE ETASQ 
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=Race Regression Range Race*Regression Race*Range Regression*Range Race*Regression*Range
  /DESIGN=Group.

DATASET ACTIVATE DataSet3.
GLM Dprime.Asian.biasm.15 Dprime.Asian.biasm.21 Dprime.Asian.biasm.27 Dprime.Asian.biasm.33 
    Dprime.Asian.biasm.39 Dprime.Asian.biasm.45 Dprime.Asian.biasp.15 Dprime.Asian.biasp.21 
    Dprime.Asian.biasp.27 Dprime.Asian.biasp.33 Dprime.Asian.biasp.39 Dprime.Asian.biasp.45 
    Dprime.Caucasian.biasm.15 Dprime.Caucasian.biasm.21 Dprime.Caucasian.biasm.27 
    Dprime.Caucasian.biasm.33 Dprime.Caucasian.biasm.39 Dprime.Caucasian.biasm.45 
    Dprime.Caucasian.biasp.15 Dprime.Caucasian.biasp.21 Dprime.Caucasian.biasp.27 
    Dprime.Caucasian.biasp.33 Dprime.Caucasian.biasp.39 Dprime.Caucasian.biasp.45
  /WSFACTOR=Race 2 Polynomial Regression 2 Polynomial Range 6 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(Range*Regression*Race Race*Regression) TYPE=BAR ERRORBAR=NO MEANREFERENCE=NO
  /EMMEANS=TABLES(OVERALL) 
  /EMMEANS=TABLES(Race) COMPARE ADJ(BONFERRONI)
  /EMMEANS=TABLES(Regression) COMPARE ADJ(BONFERRONI)
  /EMMEANS=TABLES(Range) COMPARE ADJ(BONFERRONI)
  /EMMEANS=TABLES(Race*Regression) 
  /EMMEANS=TABLES(Race*Range) 
  /EMMEANS=TABLES(Regression*Range) 
  /EMMEANS=TABLES(Race*Regression*Range) 
  /PRINT=DESCRIPTIVE ETASQ 
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=Race Regression Range Race*Regression Race*Range Regression*Range Race*Regression*Range.    
