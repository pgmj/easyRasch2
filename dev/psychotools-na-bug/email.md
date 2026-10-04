To: Achim Zeileis <Achim.Zeileis@R-project.org>
Subject: psychotools 0.7.6: pcmodel()/rsmodel() wrong estimates when no response pattern is complete

Dear Achim,

I maintain easyRasch2, which uses psychotools for CML estimation, and I think
I have found a bug in the missing-data branch of pcmodel() and rsmodel() in
psychotools 0.7.6.

**Symptom.** When no row of the response matrix is complete and every missing
data pattern has the same number of items and categories, both functions
return a positive log-likelihood, item parameters that are off by up to 0.6
logits in my example, and a warning that the Hessian could not be inverted.
A simple case is a respondent sample where everyone skipped exactly one item.
I ran into it with item-split data for DIF analysis, where an item is split
into one column per group and no respondent can have both. Adding a handful
of complete rows makes the problem disappear, which is why it is easy to
miss. raschmodel() is not affected in my tests.

**Cause.** The per-pattern parameter lists are built as

    mapply(split, esf_par, mapply(function(x, y) rep(1:x, y), m_i, oj_max_i))

in cloglik(), agrad() and the `full` ESF computation of pcmodel(), and in the
corresponding three places in rsmodel() (with rep.int and oj_i). When all
patterns have the same shape, both mapply() calls simplify their results from
lists to matrices, so elementary_symmetric_functions() is given the wrong
structure. A complete-case pattern has a different length, which keeps the
results as lists.

**Fix.** Adding SIMPLIFY = FALSE to both mapply() calls in all three places
resolves it. With the patched functions:

- the no-complete-row example gives a sensible log-likelihood and estimates
  within 0.04 logits of the complete-data fit,
- thresholds agree with eRm::PCM() on the same incomplete data to about 2e-5
  after centring, and a split-item DIF example agrees with eRm::PCM() and
  eRm::RSM() to three decimals,
- coefficients are identical to the unpatched functions on complete data and
  on data with 10% of cells missing at random.

The attached reprex.R reproduces the problem with psychotools alone (eRm is
optional, as an independent check) and applies the patch at run time by
editing the deparsed function, so the change is easy to inspect.

Thank you for psychotools. It is the CML engine behind much of easyRasch2.

Best regards,
Magnus Johansson
Karolinska Institutet

Attachment: reprex.R
