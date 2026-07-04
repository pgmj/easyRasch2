# RMitemInfitCutoffMI errors on a non-mids object

    Code
      RMitemInfitCutoffMI(data.frame(a = 0:1, b = 1:0))
    Condition
      Error:
      ! `mids_object` must be a 'mids' object as returned by mice::mice().

# RMitemInfitMI errors on a non-mids object

    Code
      RMitemInfitMI(data.frame(a = 0:1, b = 1:0))
    Condition
      Error:
      ! `mids_object` must be a 'mids' object as returned by mice::mice().

# RMitemInfitMI rejects malformed cutoff arguments

    Code
      RMitemInfitMI(imp, cutoff = "abc")
    Condition
      Error:
      ! `cutoff` must be NULL, the return value of RMitemInfitCutoff() / RMitemInfitCutoffMI(), or its $item_cutoffs data.frame.

---

    Code
      RMitemInfitMI(imp, cutoff = data.frame(Item = "I1", infit_low = 0.8))
    Condition
      Error:
      ! `cutoff` data.frame is missing required columns: infit_high.

