A6 raw-correlation version, Eric's pipeline (lib/eric_pass.py, lib/eric_xsbj.py).
x:rawmean = unweighted mean of raw per-subject Pearson r (no per-subject ranking), edges = top fraction by that value (ranked once at the end). Replicability 31.3 vs rank-averaged 37.8 at 1%; preservation 10.2 vs 11.05; at 5%: 4.5 vs 5.6.
Sign-based edges (positive = mean raw r > 0; the reviewer's "signed Pearson"): replicability of x:rawsign 1.23 (19 of 216 pairs <=1) vs xSubject (s:std, r>0) 1.50 (none <=1); preservation 1.44. About half of all gene pairs are positive, so sign agreement is weakly informative and xCell is not better by this definition.
Not possible: a rank-free version of xSubject/xBulk edges differs from the rank version only by the sign rule (ordering identical).
