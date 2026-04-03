# Repeated Measures ANOVA – SPSS Syntax Documentation

## Software
IBM SPSS Statistics for Windows, Version 24.0  

Analysis performed via:  
**Power Analysis → General Linear Model → Repeated Measures**

---

## Main Analysis: Mixed-Design ANOVA

A mixed-design ANOVA was conducted on mean **d-prime / Weibull-PSE** values.

### Factors

- **Within-subject factors:**
  - Race: Own-race / Other-race
  - Regression: Bias+ / Bias−
  - Range: 15%, 21%, 27%, 33%, 39%, 45%

- **Between-subject factor:**
  - Group: CP / TD

- **Additional analysis:**
  - Repeated for subset: participants aged ≤ 42

---

## SPSS Syntax: Full Model

```spss
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
```

---

## Follow-up Analysis: Within-Group Models

Separate analyses were conducted for each group (CP and TD independently).

### SPSS Syntax

```spss
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
```

---

## Notes

- Type III sums of squares were used (`SSTYPE(3)`).
- Polynomial contrasts were applied to all within-subject factors.
- Bonferroni correction was used for multiple comparisons.
- Effect sizes reported as partial eta squared (η²).
- Alpha level set to 0.05.

---

## Summary

This analysis tests:
- Main effects of **Race, Regression, Range, and Group**
- All interaction terms up to **4-way interaction**
- Follow-up analyses within each group to decompose interactions
